import cooling as cl

T_s = 20.0 # air temp.
T_after_death01 = 26.7 # at 2 pm
T_after_death02 = 25.8 # at 3 pm

# using FE to find temp. coeff. k :

dt = 1 # 1 hour per time step
k = (-T_after_death02 + T_after_death01)/(dt*(T_after_death02 - T_s))

print "Temperature coefficient is found to be : ",k

t_end = 250 # simulation run time
T,t = cl.cooling2(37.0, k, T_s, t_end, T_after_death02,dt, 0.5)

n=0
while T[n]>0 :
    n = n+1
# 3 pm = 15.00 hours :
print "Death occured around ",15-n-1, "am."

