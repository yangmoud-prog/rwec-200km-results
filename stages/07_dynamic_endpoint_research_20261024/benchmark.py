from calibrate import *
data=dataset('A',0);v=np.array([0.,np.log(4.),1.]);batch_loss(v,*data,0)
t=time.perf_counter()
for _ in range(20):r=batch_loss(v,*data,0)
print('objective_seconds',(time.perf_counter()-t)/20,'loss',r)
