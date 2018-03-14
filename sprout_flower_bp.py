def write_H1_Nterm(filename):
    H1_nterm='1 S .\n\
2 K .\n\
3 D .\n\
4 T .\n\
5 E .\n\
6 D .\n\
7 S .\n\
8 R .\n\
9 K .\n\
10 I .\n\
11 W .\n\
12 E .\n\
13 D .\n\
14 I .\n\
15 R .\n\
16 R .\n\
17 L .\n\
18 L .\n\
19 E .\n\
20 E .\n\
21 x H\n'
    out = open(filename, 'a')
    out.write(H1_nterm)

def write_H2_Cterm(filename):
    out = open(filename, 'a')
    C1_cterm='22 x H\n\
23 I .\n\
24 A .\n\
25 E .\n\
26 M .\n\
27 L .\n\
28 V .\n\
29 R .\n\
30 I .\n\
31 A .\n\
32 E .\n\
33 L .\n\
34 L .\n\
35 S .\n\
36 R .\n\
37 Q .\n\
38 T .\n\
39 E .\n\
40 Q .\n\
41 R .'
    out.write(C1_cterm)
def write_H_extension(filename, length):
    out = open(filename, 'a')
    helix_stub='0 x H\n'
    count = 0
    while count < length:
        out.write(helix_stub)
        count += 1

def write_loop(filename, length):
    out = open(filename,'a')
    loop_stub='0 x L\n'
    count = 0
    while count < length:
        out.write(loop_stub)
        count += 1

def write_blueprint(i,j,k):
    write_H1_Nterm('test_'+str(i)+'_'+str(j)+'_'+str(k)+'.bp')
    write_H_extension('test_'+str(i)+'_'+str(j)+'_'+str(k)+'.bp',i)
    write_loop('test_'+str(i)+'_'+str(j)+'_'+str(k)+'.bp',j)
    write_H_extension('test_'+str(i)+'_'+str(j)+'_'+str(k)+'.bp',k)
    write_H2_Cterm('test_'+str(i)+'_'+str(j)+'_'+str(k)+'.bp')

def write_all_blueprint(i,x,j,y,k,z):
    count1 = i
    while count1 <= x:
        count2  = j
        while count2 <= y:
            count3 = k
            while count3 <= z:
                print ('test_'+str(count1)+'_'+str(count2)+'_'+str(count3)+'.bp')
                write_blueprint(count1,count2,count3)
                count3 +=1
            count2 +=1
        count1 += 1

write_all_blueprint(15,20,2,5,15,20)
