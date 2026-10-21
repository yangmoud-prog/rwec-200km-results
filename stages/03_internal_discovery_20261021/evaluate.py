"""Independent full-process evaluation, including paired responses and fixed-step checks."""
import argparse,math,re
from pathlib import Path
import numpy as np
import pandas as pd
from prepare import OUT,ROOT,read,save,frame,phase
from fit import load,hy,core

def configs():
    result=[dict(id='RWEC0',label='RWEC-0',type='base'),dict(id='RWEC1',label='RWEC-1',type='empirical',**read(OUT/'fits/empirical.json'))]
    if (OUT/'fits/jointbase.json').exists():result.append(dict(id='RWEC0J',label='RWEC-0均衡率定',type='jointbase',base_fits=read(OUT/'fits/jointbase.json')))
    for mode in [1,2]:result.append(dict(id=f'P{mode}',label=f'P{mode}状态修正',type='physical',mode=mode,**read(OUT/'fits'/f'P{mode}.json')))
    gate=read(OUT/'fits/neural_gate.json')
    seed=gate['chosen_seed_by_training_loss']
    for s in [20261021,20261022]:
        result.append(dict(id=f'NN{s}',label='小网络' if s==seed else '小网络复核',type='neural',primary=s==seed,**read(OUT/'fits'/f'neural_{s}.json')))
    if (OUT/'fits/symbolic_candidates.json').exists():
        for c in read(OUT/'fits/symbolic_candidates.json')['candidates']:
            result.append(dict(label=f'解析候选{c["id"][-1]}',type='symbolic',**c))
    return result

def metric(q,y,t):
    q=np.asarray(q);y=np.asarray(y);dif=q-y;den=np.sum((y-y.mean())**2)
    nse=1-float(np.sum(dif*dif)/den) if den>1e-15 else np.nan
    if np.std(q)>1e-12 and np.std(y)>1e-12 and abs(y.mean())>1e-12:
        kge=1-math.sqrt((np.corrcoef(q,y)[0,1]-1)**2+(np.std(q)/np.std(y)-1)**2+(q.mean()/y.mean()-1)**2)
    else:kge=np.nan
    return dict(NSE=nse,KGE=float(kge),RMSE=float(np.sqrt(np.mean(dif*dif))),bias_m3s=float(dif.mean()),
        peak_error_m3s=float(q.max()-y.max()),peak_time_error_h=float(t[np.argmax(q)]-t[np.argmax(y)]),
        min_prediction_m3s=float(q.min()),negative_predictions=int((q<0).sum()))

def predict(c,q,h,ns,L,f,dt=1.):
    if c['type']=='jointbase':f=c['base_fits'][int(round(L/10000))-1]
    qr,s,ib,_=core.route(q,h,len(q),f['k'],f['n'],f['tau_h'],c.get('alpha',0) if c['type']=='empirical' else 0,0,L,dt)
    if c['type'] in ['base','jointbase','empirical']:return qr,s,ib,np.zeros(len(q)),np.zeros(len(q),dtype=int),float(qr.min())
    scale=500*L*np.maximum(h-5,0)/1e6
    if c['type']=='physical':eq=scale*c['alpha'];tb=c.get('tb_h') or 0.;mode=c['mode']
    else:
        features=np.stack([np.full(len(q),L/100000),np.maximum(h-5,0)/6,ib/6000],axis=-1)
        if c['type']=='neural':
            th=np.array(c['theta']);coef,*_=hy.network(th,features);tb=math.exp(th[21])
        else:
            import sympy as sp
            expression=sp.sympify(c['equation_logit']);fn=sp.lambdify(sp.symbols('l h i'),expression,'numpy')
            # The coefficient is multiplied by a known zero at H=H0; its
            # algebraic value outside the active-stage domain is not required.
            active_h=scale>0;g=np.zeros_like(q)
            with np.errstate(all='ignore'):g[active_h]=fn(*features[active_h].T)
            if not np.all(np.isfinite(g)):raise ValueError(f'Nonfinite symbolic expression {c["id"]}')
            coef=2/(1+np.exp(-np.clip(g,-700,700)));tb=c['tb_h']
        eq=scale*coef;mode=2
    pp,bb,aa,raw=hy.response(qr[None,:],eq[None,:],np.array([len(q)],dtype=np.int64),dt,tb,mode,True)
    # Independent, unconstrained original trajectory; never replace its negatives.
    unbounded,*_=hy.response(qr[None,:],eq[None,:],np.array([len(q)],dtype=np.int64),dt,tb,mode,False)
    return pp[0],s,ib,bb[0],aa[0],float(unbounded.min())

def full():
    d,meta=load();fits=read(OUT/'fits/base.json');models=configs();records=[];phases=[];balances=[];preds={};fail=[]
    for c in models:
        pred=np.zeros_like(d['qin']);states=pred.copy();extras=pred.copy();activations=np.zeros_like(pred,dtype=np.int64)
        for l,f in enumerate(fits):
            L=(l+1)*10000.
            for j,N in enumerate(d['ns']):
                t=np.arange(N)-48.;m=meta[j]
                try:
                    q,s,ib,b,act,unmin=predict(c,d['qin'][l,j,:N],d['h'][j,:N],N,L,f)
                    assert np.all(np.isfinite(q))
                except Exception as e:
                    fail.append(dict(model=c['id'],case=m['case'],Lr_km=L/1000,error=str(e)));continue
                pred[l,j,:N]=q;pred[l,j,N:]=q[-1];states[l,j,:N]=s;extras[l,j,:N]=b;activations[l,j,:N]=act
                row=dict(model=c['id'],case=m['case'],event_id=m['event_id'],subset=m['subset'],role=m['role'],kind=m['kind'],
                    Lr_km=int(L/1000),**metric(q,d['obs'][j,:N],t),
                    constrained_steps=int((act!=0).sum()),constrained_fraction=float((act!=0).sum()/max(1,N-1)),
                    unconstrained_original_min_m3s=unmin)
                records.append(row)
                r,phys,transit=hy.balance(d['qin'][l,j,:N],s,ib,b,q,1,f['tau_h'])
                balances.append(dict(model=c['id'],case=m['case'],event_id=m['event_id'],Lr_km=int(L/1000),
                    incremental_state_balance_max_m3=float(abs(r).max()),physical_control_balance_max_m3=float(abs(phys).max()),
                    max_abs_transit_volume_m3=float(abs(transit).max())))
                labels=np.array([phase(tt,m['event']) for tt in t])
                for label in np.unique(labels):
                    mask=labels==label;dif=q[mask]-d['obs'][j,:N][mask]
                    phases.append(dict(model=c['id'],case=m['case'],event_id=m['event_id'],subset=m['subset'],role=m['role'],kind=m['kind'],
                        Lr_km=int(L/1000),phase=label,samples=int(mask.sum()),RMSE=float(np.sqrt(np.mean(dif*dif))),bias_m3s=float(dif.mean())))
        np.savez_compressed(OUT/'predictions'/f'{c["id"]}.npz',q=pred,s=states,b=extras,active=activations)
        preds[c['id']]=pred
        print('完整过程评价',c['id'],flush=True)
    frame(OUT/'tables/case_metrics.csv',records);frame(OUT/'tables/phase_metrics.csv',phases);frame(OUT/'tables/balance.csv',balances)
    save(OUT/'fits/evaluation_failures.json',fail)
    if fail:raise RuntimeError('Candidate runtime failures retained; review before aggregate reporting.')
    M=pd.DataFrame(records)
    cols=['NSE','KGE','RMSE','bias_m3s','peak_error_m3s','peak_time_error_h','min_prediction_m3s','constrained_fraction']
    # Main-event means average paired conditions and nested lengths within flood first.
    main=M[M.role=='main'];events=main.groupby(['model','event_id','subset'])[cols].mean().reset_index()
    frame(OUT/'tables/event_metrics.csv',events)
    group=events.groupby(['model','subset'])[cols].agg(['mean','median','min','max']);group.columns=['_'.join(x) for x in group.columns]
    group.reset_index().to_csv(OUT/'tables/event_summary.csv',index=False,encoding='utf-8-sig')
    for filename,df,keys in [('condition_summary.csv',main,['model','subset','kind']),
        ('length_summary.csv',main,['model','subset','kind','Lr_km']),
        ('control_summary.csv',M[M.role!='main'],['model','subset'])]:
        g=df.groupby(keys)[cols].agg(['mean','median','min','max']);g.columns=['_'.join(x) for x in g.columns]
        g.reset_index().to_csv(OUT/'tables'/filename,index=False,encoding='utf-8-sig')
    base=events[events.model=='RWEC0'].set_index(['event_id','subset'])
    gains=[]
    for _,r in events.iterrows():
        b=base.loc[(r.event_id,r.subset)]
        gains.append(dict(model=r.model,event_id=r.event_id,subset=r.subset,NSE_gain=r.NSE-b.NSE,RMSE_change=r.RMSE-b.RMSE))
    frame(OUT/'tables/event_improvements.csv',gains)
    pair_response(d,meta,preds)
    save(OUT/'checkpoint.json',dict(stage='full_evaluation_complete',test_read=0,models=[c['id'] for c in models]))

def pair_response(d,meta,preds):
    rows=[]
    def one(a,b,category,factor,description):
        N=min(d['ns'][a],d['ns'][b]);t=np.arange(N)-48.;y=d['obs'][a,:N]-d['obs'][b,:N]
        for l in range(10):
            di=d['qin'][l,a,:N]-d['qin'][l,b,:N]
            for name,pr in preds.items():
                q=pr[l,a,:N]-pr[l,b,:N]
                rows.append(dict(model=name,event_id=meta[a]['event_id'],subset=meta[a]['subset'],category=category,
                    factor=factor,case_A=meta[a]['case'],case_B=meta[b]['case'],difference_definition=description,
                    Lr_km=(l+1)*10,common_start_h=-48.,common_end_h=float(t[-1]),
                    input_difference_RMSE_m3s=float(np.sqrt(np.mean(di*di))),
                    target_difference_min_m3s=float(y.min()),target_difference_max_m3s=float(y.max()),
                    **metric(q,y,t)))
    for eid in sorted({m['event_id'] for m in meta if m['role']=='main'}):
        a=next(j for j,m in enumerate(meta) if m['event_id']==eid and m['role']=='main' and m['kind']=='var')
        b=next(j for j,m in enumerate(meta) if m['event_id']==eid and m['role']=='main' and m['kind']=='ref')
        one(a,b,'main_pair','变动减固定','工况A（变动水位）流量减工况B（固定5m水位）流量；负值表示变动边界抑流')
    bg=sorted({m['event_id'] for m in meta if m['role']!='main'})
    for eid in bg:
        group={int(re.search(r'_C(\d+)_',m['case']).group(1)):j for j,m in enumerate(meta) if m['event_id']==eid and m['role']!='main'}
        for name,bit,desc in [('水位幅度',4,'工况A（高幅水位）流量减工况B（低幅水位）流量'),
            ('上涨时机',2,'工况A（较晚上涨）流量减工况B（较早上涨）流量'),
            ('维持时间',1,'工况A（维持0.20D）流量减工况B（不维持）流量')]:
            for z in range(8):
                if z&bit:one(group[z+1],group[(z^bit)+1],'factorial_control',name,desc)
    frame(OUT/'tables/paired_response_metrics.csv',rows)
    R=pd.DataFrame(rows);g=R.groupby(['model','subset','category','factor','Lr_km'])[['NSE','KGE','RMSE']].agg(['mean','median','min','max'])
    g.columns=['_'.join(c) for c in g.columns];g.reset_index().to_csv(OUT/'tables/paired_response_summary.csv',index=False,encoding='utf-8-sig')

def sensitivity():
    d,meta=load();fits=read(OUT/'fits/base.json');models=[c for c in configs() if c['type']!='neural' or c['primary']]
    rows=[];ex={};bgs=['E023','E006','E007','E069']
    for eid in bgs:
        j=next(j for j,m in enumerate(meta) if m['role']=='main' and m['event_id']==eid and m['kind']=='var')
        N=d['ns'][j];t=np.arange(N,dtype=float)
        for l in [0,4,9]:
            for c in models:
                paths={}
                for dt in [1.,.5,.25]:
                    tt=np.arange((N-1)*int(1/dt)+1)*dt
                    assert np.all(np.diff(t)>0)
                    qi=np.interp(tt,t,d['qin'][l,j,:N]);hh=np.interp(tt,t,d['h'][j,:N])
                    q,*_=predict(c,qi,hh,len(qi),(l+1)*10000.,fits[l],dt)
                    step=int(round(1/dt));agg=np.r_[q[0],q[1:].reshape(N-1,step).mean(axis=1)]
                    paths[dt]=agg
                    ex[f'{eid}_{l}_{c["id"]}_{str(dt)}']=agg
                for dt in [.5,.25]:
                    dif=paths[dt]-paths[1.]
                    rows.append(dict(model=c['id'],event_id=eid,subset=meta[j]['subset'],Lr_km=(l+1)*10,dt_h=dt,comparison='vs_1h',
                        RMSE_vs_1h_m3s=float(np.sqrt(np.mean(dif*dif))),max_abs_vs_1h_m3s=float(abs(dif).max()),
                        RMSE_vs_hec_m3s=float(np.sqrt(np.mean((paths[dt]-d['obs'][j,:N])**2))),
                        peak_change_m3s=float(paths[dt].max()-paths[1.].max())))
                fine=paths[.25]-paths[.5]
                rows.append(dict(model=c['id'],event_id=eid,subset=meta[j]['subset'],Lr_km=(l+1)*10,dt_h=.25,comparison='0.25h_vs_0.5h',
                    RMSE_vs_1h_m3s=float(np.sqrt(np.mean(fine*fine))),max_abs_vs_1h_m3s=float(abs(fine).max()),
                    RMSE_vs_hec_m3s=np.nan,peak_change_m3s=np.nan))
    frame(OUT/'tables/timestep_sensitivity.csv',rows)
    np.savez_compressed(OUT/'predictions/timestep_examples.npz',**ex)
    # Conditional training-only parameter profile; no refitting core or validation.
    from fit import arrays,mse
    A=arrays(d,meta,'train');prof=[]
    for tb in np.geomspace(.1,120,21):
        for alpha in np.linspace(0,.4,21):
            q,*_=hy.response(A['qr'],A['scale']*alpha,A['ns'],1.,tb,2,True)
            prof.append(dict(alpha=float(alpha),tb_h=float(tb),train_loss=mse(q,A)))
    frame(OUT/'tables/P2_parameter_profile.csv',prof)
    params=[]
    for l,f in enumerate(fits):
        params.append(dict(Lr_km=(l+1)*10,k=f['k'],n=f['n'],tau_h=f['tau_h'],
            k_lower_hit=f['k']<=1e-6*(1+1e-5),k_upper_hit=f['k']>=100*(1-1e-5),
            n_lower_hit=f['n']<=.01*(1+1e-5),n_upper_hit=f['n']>=300*(1-1e-5),
            tau_lower_hit=f['tau_h']<=.001,tau_upper_hit=f['tau_h']>=149.999,
            k_near_lower=f['x'][0]<=.01,k_near_upper=f['x'][0]>=.99,
            n_near_lower=f['x'][1]<=.01,n_near_upper=f['x'][1]>=.99,
            tau_near_lower=f['x'][2]<=.01,tau_near_upper=f['x'][2]>=.99,
            training_loss=f['loss_after'],nfev=f['total_nfev']))
    frame(OUT/'tables/base_parameters.csv',params)
    print('固定参数步长检查与训练参数剖面完成',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['full','sensitivity']);a=p.parse_args()
    {'full':full,'sensitivity':sensitivity}[a.action]()
