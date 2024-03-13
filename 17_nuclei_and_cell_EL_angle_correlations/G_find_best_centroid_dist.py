


def readfile(filename):
    d={}
    for line in filename:
        l=line.split()
        d[l[5]]=1
        d[l[6]]=1
    return d


def main():
    f1=open('cell/spherical_coordinate.txt')
    f2=open('nuc/spherical_coordinate.txt')

    dcel=readfile(f1)
    dnuc=readfile(f2)

    f=open('cell_nuc_pairs_olf.dat')
    cont=f.readlines()

    close1=[]
    close2=[]
    for j in range(len(cont)):
        l=cont[j].split()
        #print(l[0])
        if l[0] in dcel:
            close1.append(float(l[2]))
        if l[1] in dnuc:
            close2.append(float(l[2]))

    #for key in dcel:
    a=sorted(close1)
    b=sorted(close2)
    a=a[::-1]
    b=b[::-1]

    print(a[0:10])
    print(b[0:10])


main()
