
import sys, math, os, numpy

# ADD the PDB file directory HERE
PDB_DIR = '/Users/soubiaso/nmrdata/PRE_backCal/data_base_calculation/PDBs/'

R2_red  = 4.0 # not quite understand
t_INEPT = 0.01

K       = 1.23e-32 
tauC = 34*1e-9 #74*1e-9
tc      = 1.0/(1.0/(tauC) + 1.0/(100*1e-9)) # s

rh = 2.6752*1e8        
w       = -rh*18.8


#### internsity ratio #####

Delay = 3.6*2/1000.0 # in second
R2    = 42 #42  # per second

def read_r2(fn):
    with open(fn) as f:
        lines = f.readlines()
    data = {}
    for l in lines:
        l = l.split()
        if l:
            resi, r2 = int(l[0]), float(l[1])
            data[resi] = r2
    return data


def inten_r(gamma2=0.0, r2_0=R2):

    
    return (r2_0*math.exp(-gamma2*Delay))/ (r2_0+gamma2)

def cal(r=0.0, tauR=100.0):
    tc      = 1.0/(1.0/(tauC) + 1.0/(tauR*1e-9))
    r2 = (K*(4*tc + (3*tc)/(1+(w*tc)**2 ) ) )/(r**6)
    
    return r2

def dist(a1,a2):
    x1, y1, z1 = a1
    x2, y2, z2 = a2
    
    dist = math.sqrt((x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0)
    return dist

def read_pdb(fn):
    p_posi  = []
    c5_posi = []
    pro_ch3 = []
    n_helix = [[], [], []]
    
    with open(fn) as f:
        lines = f.readlines()

    for l in lines:
        ll = l.split()
        if ll:
            if ll[0] == 'ATOM' and ll[-1]=='MEMB':
                
                if l[13:16] in ['C25']:
                    x, y, z = float(l[31:38]), float(l[38:46]), float(l[46:54])
                    c5_posi.append(([x,y,z],l[23:26]))
                    

            elif ll[0] == 'ATOM' and ll[-1] == 'PROA':
                if ll[3] == 'ILE' and ll[2] == 'CD': pass
                elif ll[3] == 'LEU' and ll[2] == 'CD1': pass
                elif ll[3] == 'VAL' and ll[2] == 'CG1': pass
                elif int(l[22:26])<12 and l[13:15] == 'CA': pass
                else:
                    continue
                
                x, y, z = float(l[31:38]), float(l[38:46]), float(l[46:54])
                if ll[2] != 'CA':
                    pro_ch3.append(([x,y,z], ll[3]+'-'+ll[5]))
                else:
                    n_helix[0].append(x)
                    n_helix[1].append(y)
                    n_helix[2].append(z)
                #print [x,y,z], ll[3]+'-'+ll[5]
    center_h = sum(n_helix[0])/float(len(n_helix[0])),\
               sum(n_helix[1])/float(len(n_helix[1])),\
               sum(n_helix[2])/float(len(n_helix[2]))

    return  c5_posi, pro_ch3, center_h



def main(fn,tmp={}, r2_0={}):

    print ('# %s' % fn )
    c5, ch3, center_h = read_pdb(fn)

    pre = {}
    pre1 ={}
    tmpp = []
    new_c5 = []
    for i,j in enumerate(c5):
        tmpp.append( dist(j[0],center_h) )

    numLipid = len(tmpp)
    
    for c in sorted(range(len(tmpp[:numLipid//2])), key=lambda k: tmpp[:numLipid//2][k])[:50]:
        new_c5.append(c5[c])
            
    for c in sorted(range(len(tmpp[numLipid//2:])), key=lambda k: tmpp[numLipid//2:][k])[:50]:
        new_c5.append(c5[numLipid//2:][c])
 
    for i in new_c5:
   
        for j, k in ch3:
            d =  dist(i[0],j)
            
      
            r2 = cal(d*1e-8)
   
            if r2 >100:
                r2 = 100
    
            k = int(k[4:])
            
            if k not in pre.keys():
                pre[k] = []
                pre1[k] = []
            pre[k].append(cal(d/1e8))
            pre1[k].append(cal(d/1e8, 500))
            
        
    for x in  sorted(pre.keys()):
        avg = sum(pre[x])/ float(len(pre[x]))
        avg1 = sum(pre1[x])/ float(len(pre1[x]))

        if x not in r2_0.keys():
            r2_0[x] = 41.69346613

        iR = inten_r(avg, r2_0[x])
        iR1 = inten_r(avg1, r2_0[x])


        if x not in tmp.keys():
            tmp[x] = []
        tmp[x].append( inten_r(avg, r2_0[x]))
        print (x, inten_r(avg), inten_r(avg1),inten_r(avg, r2_0[x]),  iR1 )
    print ()


    return tmp

if __name__ == '__main__':

    pdbDir = [1]
    try:
        r2_0 = read_r2('./r2.dat')
        files = os.listdir(PDB_DIR)
    except:
        print ("""
########################  HOW TO USE  ########################
#
#    python PRE-db-calculation.py
#
# Two files need to set up before run:   
#    1) r2.dat   The file name needs exact the same. 
#       Check example r2.dat data format.
#    2) PDB_DIR variable in line #5 needs to provide.
#
########################     END      ########################
""")

    data  = {}

    for ii in pdbDir[:]:  
        files = os.listdir(PDB_DIR)
        data  = {}
        for fn in files[:]:
            if '.pdb' in fn:
                data = main(PDB_DIR+fn,data)
