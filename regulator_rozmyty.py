import matplotlib.pyplot as p
import pylab
import numpy as np
import skfuzzy as fuzz

class zbiornik:
    """docstring for [object Object]."""
    def __init__(self,hmax,h,s,s0):
        self.hmax = hmax
        self.h = h
        self.s = s
        self.s0 = s0
        self.vmax = self.s*self.hmax
        self.v = self.h*self.s
    def a_h(self): #akutalny poziom wody h
        self.h = self.v/self.s
        return self.h
    def aktualna_obietosc(self,Q1,Q2,Q3):
        self.v += Q1+Q2-Q3
        return self.v

class PID:
    def __init__(self, kp,Ti,Td):
        self.kp = kp
        self.Ti = Ti
        self.Td = Td
def draw(nazwa,y,T,i):
    pylab.plot(T, y, 'C'+str(i), label=nazwa)
    pylab.xlabel('Czas próbkowania')
    pylab.ylabel('Wyjście')
    pylab.title('Obiekt inercyjny')
    pylab.legend()
    pylab.draw()

def draw_save():
    pylab.grid(which = 'both',alpha=1)
    pylab.show()
    #fig1.savefig('img.png', dpi=100)
def regulacja_przyrostowa(ob,r):
    suma_e = 0.0
    for k in range(1, 101, 1):
        y_ca.append(Q1[k]/ob.v*(Ca1[k]-y_ca[k])+Q2[k]/ob.v*(Ca2[k]-y_ca[k]))
        w.append(w[k-1])
        e_temp = w[k-1]-y_ca[k-1]
        e_temp = w[k-1]-y[k-1]
        e.append(e_temp)
        suma_e += e[k]
        Q1.append(r.kp*(e[k] + Tp/r.Ti*suma_e + r.Td/Tp*(e[k] - e[k-1])))

y_ca = [0.75, ] # uchyb
w = [2.0, ] #wartość zadana Ca
e = [w[0] - y_ca[0], ] #wyjscie
suma_e = 0.0
q1 = [0.0, ] # sterujemy


T = [0, ] #czas
Tp = 1
ca1 = [0.2, ]
ca2 = [0.6, ]
q2 = [1.0, ] # wartość stała zakłócenie stałe
alfa = 0.3
g = 9.81
ob = zbiornik(3,3,1,0.1)
Q3 = alfa*ob.s0*np.sqrt(2*g*ob.a_h())

r = PID(0.8,0.09,0.1)
#zakresy wartosci uchybu, zmiany uchybu, wielkosci sterujacej
x_e = np.arange(-1.5, 2, 0.5)
x_ce = np.arange(-3, 4, 1)
x_u  = np.arange(-4, 5, 1)
#funkcje przynaleznosci zbiorow rozmytych
e_du = fuzz.trimf(x_e, [-1.5, -1.5, -1])
e_su = fuzz.trimf(x_e, [-1.5, -1, -0.5])
e_mu = fuzz.trimf(x_e, [-1, -0.5, 0])
e_z = fuzz.trimf(x_e, [-0.5, 0, 0.5])
e_md = fuzz.trimf(x_e, [0, 0.5, 1])
e_sd = fuzz.trimf(x_e, [0.5, 1, 1.5])
e_dd = fuzz.trimf(x_e, [1, 1.5, 1.5])

ce_du = fuzz.trimf(x_ce, [-3, -3, -2])
ce_su = fuzz.trimf(x_ce, [-3, -2, -1])
ce_mu = fuzz.trimf(x_ce, [-2, -1, 0])
ce_z = fuzz.trimf(x_ce, [-1, 0, 1])
ce_md = fuzz.trimf(x_ce, [0, 1, 2])
ce_sd = fuzz.trimf(x_ce, [1, 2, 3])
ce_dd = fuzz.trimf(x_ce, [2, 3, 3])

u_bdu = fuzz.trimf(x_u, [-4, -4, -3])
u_du = fuzz.trimf(x_u, [-4, -3, -2])
u_su = fuzz.trimf(x_u, [-3, -2, -1])
u_mu = fuzz.trimf(x_u, [-2, -1, 0])
u_z = fuzz.trimf(x_u, [-1, 0, 1])
u_md = fuzz.trimf(x_u, [0, 1, 2])
u_sd = fuzz.trimf(x_u, [1, 2, 3])
u_dd = fuzz.trimf(x_u, [2, 3, 4])
u_bdd = fuzz.trimf(x_u, [3, 4, 4])


for n in range(0, 201, 1): #czas
    T.append(n)

for k in range(0, 201, 1):
    q3 = alfa*ob.s0*np.sqrt(2*g*ob.a_h())
    ob.aktualna_obietosc(q1[k],q2[k],q3)
    y_ca.append(((q1[k]/ob.v*ca1[k]+q2[k]/ob.v*ca2[k]-y_ca[k-1])/(1+q1[k]/ob.v+q2[k]/ob.v))/1.5)
    w.append(w[k-1])
    e_temp = w[k-1]-y_ca[k-1]
    #suma_e += e[k]
    #q1.append(r.kp*(e[k] + Tp/r.Ti*suma_e + r.Td/Tp*(e[k] - e[k-1])))
    
    if e_temp/2 > 0.75: e_temp = 1.5
    if e_temp/2 < -0.75: e_temp = -1.5
    #baza regul
    e_level_du = fuzz.interp_membership(x_e, e_du, e_temp)
    e_level_su = fuzz.interp_membership(x_e, e_su, e_temp)
    e_level_mu = fuzz.interp_membership(x_e, e_mu, e_temp)
    e_level_z = fuzz.interp_membership(x_e, e_z, e_temp)
    e_level_md = fuzz.interp_membership(x_e, e_md, e_temp)
    e_level_sd = fuzz.interp_membership(x_e, e_sd, e_temp)
    e_level_dd = fuzz.interp_membership(x_e, e_dd, e_temp)
    
    ce_level_du = fuzz.interp_membership(x_ce, ce_du, e[-1] - e_temp)
    ce_level_su = fuzz.interp_membership(x_ce, ce_su, e[-1] - e_temp)
    ce_level_mu = fuzz.interp_membership(x_ce, ce_mu, e[-1] - e_temp)
    ce_level_z = fuzz.interp_membership(x_ce, ce_z, e[-1] - e_temp)
    ce_level_md = fuzz.interp_membership(x_ce, ce_md, e[-1] - e_temp)
    ce_level_sd = fuzz.interp_membership(x_ce, ce_sd, e[-1] - e_temp)
    ce_level_dd = fuzz.interp_membership(x_ce, ce_dd, e[-1] - e_temp)
    #wyznaczanie stopni prawdziwosci przeslanek aktywnych regul sterowania;wyznaczanie stopnii prawdziwosci konkluzji aktywnych regul sterowania
    u_activation_bdu = np.amax([np.fmin(e_level_du, ce_level_du),np.fmin(e_level_du, ce_level_su),np.fmin(e_level_du, ce_level_mu),np.fmin(e_level_su, ce_level_du),np.fmin(e_level_su, ce_level_su),np.fmin(e_level_mu, ce_level_du)])
    u_rule_bdu = np.fmin(u_activation_bdu, u_bdu)
    u_activation_du = np.amax([np.fmin(e_level_z, ce_level_du),np.fmin(e_level_du, ce_level_z),np.fmin(e_level_su, ce_level_mu),np.fmin(e_level_mu, ce_level_su)])
    u_rule_du = np.fmin(u_activation_du, u_du)
    u_activation_su = np.amax([np.fmin(e_level_md, ce_level_du),np.fmin(e_level_du, ce_level_md),np.fmin(e_level_su, ce_level_z),np.fmin(e_level_z, ce_level_su),np.fmin(e_level_mu, ce_level_mu)])
    u_rule_su = np.fmin(u_activation_su, u_su)
    u_activation_mu = np.amax([np.fmin(e_level_sd, ce_level_du),np.fmin(e_level_du, ce_level_sd),np.fmin(e_level_su, ce_level_md),np.fmin(e_level_md, ce_level_su),np.fmin(e_level_mu, ce_level_z),np.fmin(e_level_z, ce_level_mu)])
    u_rule_mu = np.fmin(u_activation_mu, u_mu)
    u_activation_z = np.amax([np.fmin(e_level_dd, ce_level_du),np.fmin(e_level_du, ce_level_dd),np.fmin(e_level_su, ce_level_sd),np.fmin(e_level_sd, ce_level_su),np.fmin(e_level_mu, ce_level_md),np.fmin(e_level_md, ce_level_mu),np.fmin(e_level_z, ce_level_z)])
    u_rule_z = np.fmin(u_activation_z, u_z)
    u_activation_md = np.amax([np.fmin(e_level_dd, ce_level_su),np.fmin(e_level_su, ce_level_dd),np.fmin(e_level_mu, ce_level_sd),np.fmin(e_level_sd, ce_level_mu),np.fmin(e_level_z, ce_level_md),np.fmin(e_level_md, ce_level_z)])
    u_rule_md = np.fmin(u_activation_md, u_md)
    u_activation_sd = np.amax([np.fmin(e_level_dd, ce_level_mu),np.fmin(e_level_mu, ce_level_dd),np.fmin(e_level_z, ce_level_sd),np.fmin(e_level_sd, ce_level_z),np.fmin(e_level_md, ce_level_md)])
    u_rule_sd = np.fmin(u_activation_sd, u_sd)
    u_activation_dd = np.amax([np.fmin(e_level_dd, ce_level_z),np.fmin(e_level_z, ce_level_dd),np.fmin(e_level_md, ce_level_sd),np.fmin(e_level_sd, ce_level_md)])
    u_rule_dd = np.fmin(u_activation_dd, u_dd)
    u_activation_bdd = np.amax([np.fmin(e_level_dd, ce_level_dd),np.fmin(e_level_sd, ce_level_dd),np.fmin(e_level_md, ce_level_dd),np.fmin(e_level_dd, ce_level_sd),np.fmin(e_level_sd, ce_level_sd),np.fmin(e_level_dd, ce_level_md)])
    u_rule_bdd = np.fmin(u_activation_bdd, u_bdd)
    
  
    #print(' ')
    #agregacja stopni prawdziwosci konkluzji aktywnych regul sterowania
    aggregated = np.fmax(u_rule_bdu,np.fmax(u_rule_du, np.fmax(u_rule_su,np.fmax(u_rule_mu,np.fmax(u_rule_z,np.fmax(u_rule_md,np.fmax(u_rule_sd,np.fmax(u_rule_dd,u_rule_bdd))))))))
    
    u = fuzz.defuzz(x_u, aggregated, 'centroid')
    
    q1.append(u)
    
    e.append(e_temp)
    
    q2.append(q2[-1])
    ca1.append(ca1[-1])
    ca2.append(ca2[-1])

print(q1)
draw("uchyb",y_ca,T,0)
draw("wyjscie",e,T,1)
draw("zadana",w,T,2)

fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(8, 9))
ax0.plot(x_e, e_du, 'b', linewidth=1.5,)
ax0.plot(x_e, e_su, 'g', linewidth=1.5,)
ax0.plot(x_e, e_mu, 'r', linewidth=1.5, )
ax0.plot(x_e, e_z, 'b', linewidth=1.5, )
ax0.plot(x_e, e_md, 'g', linewidth=1.5, )
ax0.plot(x_e, e_sd, 'r', linewidth=1.5,)
ax0.plot(x_e, e_dd, 'r', linewidth=1.5, )
ax0.set_title('uchyb')
ax0.legend()


ax1.plot(x_ce, ce_du, 'b', linewidth=1.5, )
ax1.plot(x_ce, ce_su, 'g', linewidth=1.5, )
ax1.plot(x_ce, ce_mu, 'r', linewidth=1.5, )
ax1.plot(x_ce, ce_z, 'b', linewidth=1.5, )
ax1.plot(x_ce, ce_md, 'g', linewidth=1.5, )
ax1.plot(x_ce, ce_sd, 'r', linewidth=1.5, )
ax1.plot(x_ce, ce_dd, 'r', linewidth=1.5, )
ax1.set_title('zmiana uchybu')
ax1.legend()


ax2.plot(x_u, u_bdu, 'b', linewidth=1.5, )
ax2.plot(x_u, u_du, 'g', linewidth=1.5, )
ax2.plot(x_u, u_su, 'r', linewidth=1.5, )
ax2.plot(x_u, u_mu, 'b', linewidth=1.5, )
ax2.plot(x_u, u_z, 'g', linewidth=1.5, )
ax2.plot(x_u, u_md, 'r', linewidth=1.5, )
ax2.plot(x_u, u_sd, 'r', linewidth=1.5, )
ax2.plot(x_u, u_dd, 'g', linewidth=1.5, )
ax2.plot(x_u, u_bdd, 'r', linewidth=1.5, )

ax2.set_title('wartosc sterująca')
ax2.legend()

draw_save()
# for t in T:
#     Q3 = alfa*ob.s0*np.sqrt(2*g*ob.a_h())
#     print(ob.aktualna_obietosc(Q1[t],Q2[t],Q3))
#     Q2.append(Q2[-1])
#     Q1.append(Q1[-1])
#
# print(Q3)
