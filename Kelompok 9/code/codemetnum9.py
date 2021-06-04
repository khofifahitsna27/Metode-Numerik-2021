#!/usr/bin/env python
# coding: utf-8

# In[11]:


#TUGAS AKHIR METODE NUMERIK 2021
#MODUL 2 Akar-Akar Persamaan
#Oleh: KELOMPOK 9 / OSEANOGRAFI

"""
Akar-akar persamaan: Bentuk persamaan-persamaan polinomial sangat sulit atau tidak mungkin diselesaikan secara manual,
Hal ini dikarenakan penyelesaian numerik dilakukan dengan perkiraan yang berurutan (iterasi),
sehingga setiap hasil hasilnya akan lebih teliti dari perkiraan sebelumnya, yaitu dengan melakukan sejumlah prosedur iterasi
yang dianggap cukup hingga didapatkan hasil perkiraan yang mendekati hasil eksak dengan toleransi kesalahan yang diijinkan.
"""

#~~~METODE SETENGAH INTERVAL~~~
"""
Metode Setengah Interval merupakan metode analisis numerik paling sederhana diantara metode-metode analisis lainnya.
Metode ini metode iteratif yang digunakan untuk mencari akar persamaan (linier maupun non-linear) yang memiliki akar real.
"""
def Setengah_Interval(X1,X2,a,b,c,d):
    X1 = X1
    X2 = X2
    a = a
    b = b
    c = c
    d = d
    error = 1
    iterasi = 0
    while(error > 0.0001):
        iterasi +=1
        FXi = (float(a*((X1)**3)))+(float(b*((X1)**2)))+(c*X1)+d
        FXii = (float(a*((X2)**3)))+(float(b*((X2)**2)))+(c*X2)+d
        Xt = (X1+X2)/2
        FXt = (float(a*((Xt)**3)))+(float(b*((Xt)**2)))+(c*Xt)+d
        if FXi*FXt > 0:
            X1 = Xt
        elif FXi*FXt < 0:
            X2 = Xt
        else:
            print("Akar Penyelesaian: ", Xt)
        if FXt < 0:
            error = FXt*(-1)
        else:
            error = FXt
        if iterasi > 100:
            print("Angka Tak Hingga")
            break
        print(iterasi, "|", FXi, "|", FXii, "|", Xt, "|", FXt, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", Xt)
    print("Toleransi Error: ", error)

#~~~METODE INTERPOLASI LINIER~~~#
"""
Metode Interpolasi Linier hampir mirip dengan metode setengah interval. Metode ini dikenal juga dengan 
metode false position dan ada untuk menutupi kekurangan pada metode setengah interval. Metode ini 
didasarkan pada interpolasi antara dua nilai dari fungsi yang mempunyai tanda berlawanan.
"""
def Interpolasi_Linier(X1,a,b,c,d):
    X1 = X1
    X2 = X1+1
    a = a
    b = b
    c = c
    d = d
    error = 1
    iterasi = 0
    while(error > 0.0001):
        iterasi +=1
        FX1 = (float(a*((X1)**3)))+(float(b*((X1)**2)))+(c*X1)+d
        FX2 = (float(a*((X2)**3)))+(float(b*((X2)**2)))+(c*X2)+d
        Xt = X2-((FX2/(FX2-FX1)))*(X2-X1)
        FXt = (float(a*((Xt)**3)))+(float(b*((Xt)**2)))+(c*Xt)+d
        if FXt*FX1 > 0:
            X2 = Xt
            FX2 = FXt
        else:
            X1 = Xt
            FX1 = FXt
        if FXt < 0:
            error = FXt*(-1)
        else:
            error = FXt
        if iterasi > 500:
            print("Angka Tak Hingga")
            break
        print(iterasi, "|", FX1, "|", FX2, "|", Xt, "|", FXt, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", Xt)
    print("Toleransi Error: ", error)

#~~~METODE NEWTON-RHAPSON~~~
"""
Metode ini paling banyak digunakan dalam mencari akar-akar persamaan, jika perkiraan awal dari akar adalah xi,
maka suatu garis singgung dapat dibuat dari titik (xi, f (xi)).
Titik dari garis singgung tersebut memotong sumbu-x, biasanya memberikan perkiraan yang lebih dekat dari nilai akar.
Syarat dari metode ini ,syarat yang harus dipenuhi adalah bahwa taksiran awal yang diberikan harus
sedekat mungkin dengan harga eksaknya.
"""
def Newton_Rhapson(X1,a,b,c,d):
    X1 = X1
    a = a
    b = b
    c = c
    d = d
    iterasi = 0
    akar = 1
    while (akar > 0.0001):
        iterasi += 1
        Fxn = (float(a*((X1)**3)))+(float(b*((X1)**2)))+(c*X1)+d
        Fxxn = (float((a*3)*X1)**2)+(float((b*2)*X1))+(c)
        xnp1 = X1-(Fxn/Fxxn)
        fxnp1 = (a*(xnp1**3))+(b*(xnp1**2))-(c*xnp1)+d
        Ea = ((xnp1-X1)/xnp1)*100
        if Ea < 0.0001:
            X1 = xnp1
            akar = Ea*(-1)
        else:
            akar = xnp1
            print("Nilai akar adalah: ", akar)
            print("Nilai error adalah: ", Ea)
        if iterasi > 100:
            break   
        print(iterasi, "|", X1, "|", xnp1, "|", akar, "|", Ea)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", xnp1)
    print("Toleransi Error: ", akar)

#~~~METODE SECANT~~~
"""
Pada dasarnya metode ini sama dengan metode Newton-Raphson, perbedaannya hanya terletak 
pada pendekatan untuk turunan pertama dari f saja. Pendekatan f' pada metode Secant didekati dengan 
ungkapan beda hingga yang didasarkan pada taksiran akar sebelumnya (beda mundur).
"""
def Secant(X1,a,b,c,d):
    X1 = X1
    X2 = X1-1
    a = a
    b = b
    c = c
    d = d
    error = 1
    iterasi = 0
    while(error > 0.0001):
        iterasi +=1
        FX1 = (float(a*((X1)**3)))+(float(b*((X1)**2)))+(c*X1)+d
        FXmin = (float(a*((X2)**3)))+(float(b*((X2)**2)))+(c*X2)+d
        X3 = X1-((FX1)*(X1-(X2)))/((FX1)-(FXmin))
        FXplus = (float(a*((X3)**3)))+(float(b*((X3)**2)))+(c*X3)+d
        if FXplus < 0:
            error = FXplus*(-1)
        else:
            error = FXplus
        if error > 0.0001:
            X2 = X1
            X1 = X3
        else:
            print("Selesai")
        if iterasi > 500:
            print("Angka Tak Hingga")
            break
        print(iterasi, "|", FX1, "|", FXmin, "|", X3, "|", FXplus, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", X3)
    print("Toleransi Error: ", error)

#~~~METODE ITERASI~~~
"""
Metode tidak langsung atau iterasi, merupakan metode yang berbasiskan terhadap aplikasi dari langkah–langkah/algoritma
sederhana yang diulang–ulang pada sistem persamaan tersebut hingga sistem persamaan mencapai
keadaan konvergen yang merepresentasikan solusi dari sistem persamaan tersebut.
"""
def Iterasi(X1,a,b,c,d):
    X1 = X1
    a = a
    b = b
    c = c
    d = d
    error = 1
    iterasi = 0
    while (error > 0.0001):
        iterasi +=1
        Fxn = (float(a*((X1)**3)))+(float(b*((X1)**2)))+(c*X1)+d
        X2 = (((float(-b*(X1**2)))-(float(c*X1))-(d))/a)**(1/3)
        Ea = (((X2-X1)/(X2))*100)
        if Ea < error:
            X1 = X2
            if Ea > 0:
                error = Ea
            else:
                error = Ea*(-1)
        else:
            error = Ea
        if iterasi > 100:
            print("Angka Tak Hingga")
            break
        print(iterasi, "|", X1, "|", X2, "|", Ea, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", X2)
    print("Toleransi Error: ", error)

#MODUL 3 Sistem Persamaan Linier dan Matriks
#Oleh: KELOMPOK 9 / OSEANOGRAFI

"""
Sistem persamaan linier (SPL) merupakan sistem operasi matematis yang terdiri atas dua atau lebih persamaan linier.
Dalam sistem persamaan linear terdapat metode langsung dengan metode tidak langsung (iterasi).
Matriks adalah susunan angka–angka (sering disebut elemen–elemen) yang diatur menurut baris 
dan kolom, berbentuk persegi panjang atau persegi dan ditulis diantara dua tanda kurung yaitu ( ) atau [ ].
"""

#~~~METODE ELEMINASI GAUSS~~~
"""
Eliminasi Gauss adalah suatu cara mengoperasikan nilai-nilai di dalam matriks
menjadi matriks yang lebih sederhana dan banyak digunakan dalam penyelesaian sistem persamaan linier.
 Metode ini mengubah persamaan linear tersebut ke dalam matriks augmentasi dan mengoperasikannya.
"""
import numpy as np
def Gauss(A, B):
    print("Nilai Matriks = \n", AB, "\n")
    n = len(B)
    for i in range(n):
            A = AB[i]
            for j in range(i + 1, n):
                B = AB[j]
                m = A[i] / B[i]
                AB[j] = A - m * B

    for i in range(n - 1, -1, -1):
            AB[i] = AB[i] / AB[i, i]
            A = AB[i]
    for j in range(i - 1, -1, -1):
            B = AB[j] 
            m = A[i]/B[i]
            AB[j] = A-m*B
    x = AB[:, 3]
    print("Hasil Matriks = \n", AB, "\n")
    print("Hasil Akhir Matriks = \n", x)


#~~~METODE GAUSS JORDAN~~~
"""
Eliminasi Gauss-Jordan adalah pengembangan dari eliminasi Gauss yang hasilnya lebih sederhana lagi.
Caranya adalah dengan meneruskan operasi baris dari eliminasi Gauss sehingga menghasilkan matriks yang Eselon-baris.
Metode ini digunakan untuk mencari invers dari sebuah matriks.
"""
import numpy as np
import sys
def GaussJordan(a,n):
    #Step 1 ===> Looping untuk pengolahan metode Gauss Jordan
    print("===============Mulai Iterasi===============")
    for i in range(n):
        if a[i][i]==0:
         sys.exit("Dibagi dengan angka not (Proses tidak dapat dilanjutkan)")
        for j in range(n):
         if i !=j:
            ratio = a[j][i]/a[i][i]
            #print ('posisi nol di:[',j,i,']', nilai rasio:',ratio)
            
            for k in range (n+1):
             a[j, k]=a[j][k]-ratio*a[i][k]
            print(a)
            print("==============")
    
    #STEP 2 =====> MEMBUAT SEMUA VARIABEL (x1, x2, x3, x4,...)==1
    ax = np.zeros((n, n+1))
    for i in range(n):
        for j in range(n+1):
            ax[i, j]=a[i][j]/a[i][i]
    print("======= Akhir Iterasi =======")
    return ax

#~~~METODE GAUSS SIEDEL~~~
"""
Metode iterasi Gauss-Seidel adalah metode yang menggunakan proses iterasi hingga diperoleh nilai-nilai yang berubah-ubah dan akhirnya
relatif konstan. Metode iterasi GaussSeidel dikembangkan dari gagasan metode iterasi pada solusi persamaan tak linier. 
Metode eliminasi gauss-seidel digunakan untuk menyelesaikan SPL yang berukuran kecil.
Metode ini dapat menghasilkan jumlah iterasi lebih sedikit dibandingkan metode Jacobi.
"""
def Gauss_Siedel(a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3,r):
    a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3 = a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3
    r = r
    def x1(x2,x3):
        return(D1-(b1*x2)-(c1*x3))/a1
    def x2(x1,x3):
        return(D2-(a2*x1)-(c2*x3))/b2
    def x3(x1,x2):
        return(D3-(a3*x1)-(b3*x2))/c3
    def error(n,o):
        return((n-o)/n)*100
    ax1,ax2,ax3= 0,0,0
    tabel="{0:1}|{1:7}|{2:7}|{3:7}|{4:7}|{5:7}|{6:7}"
    print(tabel.format("i", "x1", "x2", "x3", "e1", "e2", "e3"))
    for i in range(0,r):
        if i == 0:
            print(tabel.format(i, ax1, ax2, ax3, "-", "-", "-"))
            cx1=ax1
            cx2=ax2
            cx3=ax3
        else:
            cx1=eval("{0:.3f}".format(x1(ax2,ax3)))
            cx2=eval("{0:.3f}".format(x2(cx1,ax3)))
            cx3=eval("{0:.3f}".format(x3(cx1,cx2)))
            print(tabel.format(i, cx1, cx2, cx3, "{0:.2f}".format(error(cx1, ax1)), "{0:.2f}".format(error(cx2, ax2)), "{0:.2f}".format(error(cx3, ax3))))
        ax1=cx1
        ax2=cx2
        ax3=cx3

#~~~METODE JACOBI~~~
"""
Metode iterasi Jacobi digunakan untuk menyelesaikan persamaan linier yang proporsi koefisien nol-nya besar.
Keuntungan metode ini adalah langkah penyelesaiannya yang sederhana, sedangkan kelemahannya adalah proses
iterasinya lambat, terutama untuk persamaan linear serentak dengan ordo tinggi.
"""
def Jacobi(a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3,r):
    a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3 = a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3
    r = r
    def x1(x2,x3):
        return(D1-(b1*x2)-(c1*x3))/a1
    def x2(x1,x3):
        return(D2-(a2*x1)-(c2*x3))/b2
    def x3(x1,x2):
        return(D3-(a3*x1)-(b3*x2))/c3
    def error(n,o):
        return((n-o)/n)*100
    bx1,bx2,bx3= 0,0,0
    tabel="{0:1}|{1:7}|{2:7}|{3:7}|{4:7}|{5:7}|{6:7}"
    print(tabel.format("i", "x1", "x2", "x3", "e1", "e2", "e3"))
    for i in range(0,r):
        if i == 0:
            print(tabel.format(i, bx1, bx2, bx3, "-", "-", "-"))
            cx1=bx1
            cx2=bx2
            cx3=bx3
        else:
            cx1=eval("{0:.3f}".format(x1(bx2,bx3)))
            cx2=eval("{0:.3f}".format(x2(bx1,bx3)))
            cx3=eval("{0:.3f}".format(x3(bx1,bx2)))
            print(tabel.format(i, cx1, cx2, cx3, "{0:.2f}".format(error(cx1, bx1)), "{0:.2f}".format(error(cx2, bx2)), "{0:.2f}".format(error(cx3, bx3))))
        bx1=cx1
        bx2=cx2
        bx3=cx3


#MODUL 4 Integrasi Numerik
#Oleh: KELOMPOK 9 / OSEANOGRAFI

"""
Integrasi numerik merupakan cara alternatif untuk mengintegrasikan suatu persamaan, disamping integrasi analitis.
Integrasi analitis terkadang merupakan cara integrasi yang sulit, khususnya pada persamaan – persamaan
yang kompleks dan rumit. Disamping itu, juga fungsi-fungsi yang diintegralkan tidak berbentuk analitis
melainkan berupa titik-titik data. 
"""

import numpy as np
import matplotlib.pyplot as plt

#~~~METODE TRAPESIUM SATU PIAS~~~
"""
Metode trapesium merupakan metode pendekatan integral numerik dengan persamaan polinomial order satu. Dalam metode ini
kurva lengkung dari fungsi f (x) digantikan oleh garis lurus. Pendekatan dilakukan dengan satu pias (trapesium).
"""
def Trapesium_SatuPias(A,B,a,b,c):
    import numpy as np
    import matplotlib.pyplot as plt
    A = A
    B = B
    a = a
    b = b
    c = c
    x = np.linspace(-10,10,100)
    y = a*(x**3)+b*(x**2)+c
    x1 = A
    x2 = B
    x3 = B/2
    fx1 = a*(x1**3)+b*(x1**2)+c
    fx2 = a*(x2**2)+b*(x2**2)+c
    fx3 = 5
    plt.plot(x,y)
    plt.fill_between([x1,x2],[fx1,fx2])
    plt.xlim([-20,20]); plt.ylim([-20,20]);
    plt.title('Trapesium 1 pias')
    plt.savefig('image\Trapesium_satu_pias.png')
    L = 0.5*(fx2 + fx1)*(x2 - x1)
    print("luas dengan metode trapesium 1 pias:", L)

#~~~METODE TRAPESIUM BANYAK PIAS~~~
"""
Metode trapesium banyak pias diberikan untuk menutup kekurangan metode trapesium satu pias.
Diketahui bahwa pendekatan dengan menggunakan satu pias (trapesium) menimbulkan error yang besar.
Untuk mengurangi kesalahan yang terjadi maka kurve lengkung didekati oleh sejumlah garis lurus,
sehingga terbentuk banyak pias. Semakin kecil pias yang digunakan, hasil yang didapat menjadi semakin teliti .
"""
def Trapesium_BanyakPias(A,B,N,a,b,c):
    import numpy as np
    import matplotlib.pyplot as plt
    def trapesium(f,A,B,N):
        x = np.linspace(A,B,N+1)
        y = f(x)
        y_right = y[1:] 
        y_left = y[:-1] 
        dx = (B-A)/N
        T = (dx/2)*np.sum(y_right + y_left)
        return T
    f = lambda x : ((a*(x**3))+(b*(x**2))+c)
    A = A
    B = B
    N = N
    a = a
    b = b
    c = c
    x = np.linspace(A,B,N+1)
    y = f(x)
    X = np.linspace(A,B+1,N)
    Y = f(X)
    plt.plot(X,Y)
    for i in range(N):
        xs = [x[i],x[i],x[i+1],x[i+1]]
        ys = [0,f(x[i]),f(x[i+1]),0]
        plt.fill(xs,ys,'b',edgecolor='b',alpha=0.2)
    plt.title('Trapesium banyak pias, N = {}'.format(N))
    plt.savefig('image\Trapesium_banyak_pias')
    L = trapesium(f,A,B,N)
    print(L)

#~~~METODE SIMPSON 1/3~~~
"""
Metode simpson 1/3 adalah metode yang mencocokkan polinomial derajat 2 pada tiga titik data diskrit
yang mempunyai jarak yang sama. Hampiran nilai integrasi yang lebih baik dapat ditingkatkan
dengan menggunakan polinominterpolasi berderajat yang lebih tinggi.
"""
def Simpson_1per3(A,B,a,b,c):
    A = A
    B = B
    a = a
    b = b
    c = c
    def f(x):
        return a*(x**3) + b*(x**2) + c
    def simpson1per3(x0,xn,n):
        h = (xn - x0) / n
        integral = f(x0) + f(xn)
        for i in range(1,n):
            k = x0 + i*h
            if i%2 == 0:
                integral = integral + 2 * f(k)
            else:
                integral = integral + 4 * f(k)
        integral = integral * h/3
        return integral
    hasil = simpson1per3(A, B, 2)
    print("nilai integral metode Simpson 1/3:",hasil )

#~~~METODE SIMPSON 3/8~~~
"""
Metode Simpson 3/8 diturunkan dengan menggunakan persamaan polinomial order tiga yang melalui empat titik.
"""
def Simpson_3per8(A,B,a,b,c):
    A = A
    B = B
    a = a
    b = b
    c = c
    def simpson():
        import math
    def f(x):
        return a*(x**3)+b*(x**2)+c
    def simpson1per3(x0,xn,n):
        h = (xn-x0)/n
        integral = f(x0) + f(xn)
        for i in range(1,n):
            k = x0 + 1*h
        if i%2 ==0 :
            integral = integral + 2 * f(k)
        else:
            integral = integral + 4 * f(k)
            integral = integral * 3 * h/8
            return integral 
    hasil = simpson1per3(A, B, 2)
    print("nilai integral metode simpson 3/8:", hasil)

#MODUL 5 Persamaan Diferensial Biasa
#Oleh: KELOMPOK 9 / OSEANOGRAFI

"""
suatu persamaan differensial dikategorikan bedasarkan variabel bebasnya (Independent Variable). Pada persamaan differensial 
biasa atau sering disebut juga dengan Ordinary Differential Equations (ODE) adalah 
persamaan differensial yang hanya memiliki satu variabel bebas.
"""

import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython

plt.style.use('seaborn-poster')
ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')

#~~~METODE EULER~~~
"""
Metode Euler adalah salah satu dari metode satu langkah yang paling sederhana.
Di banding dengan beberapa metode lainnya, metode ini paling kurang teliti. Metode Euler dapat 
diturunkan dari Deret Taylor. Metode ini pada dasarnya adalah merepresentasikan solusinya dengan beberapa suku deret Taylor.
"""
def Euler(h,x0,xn,y0,a,b,c,d):
    h = h
    x0 = x0
    xn = xn
    x = np.arange(x0, xn + h, h)
    y0 = y0
    a = a
    b = b
    c = c
    d = d
    G = a*(x**3) + b*(x**2) + c*x + d
    f = lambda x, y: a*(x**3) + b*(x**2) + c*x + d
    y = np.zeros(len(x))
    y[0] = y0
    
    for i in range(0, len(x) - 1):
        y[i + 1] = y[i] + h*f(x[i], y[i])

    Galat = G-y
    print("galat yang diperoleh dari metode Euler adalah:", Galat)

    Judul = ("Grafik Pendekatan Persamaan Differensial Biasa Dengan Metode Euler")
    plt.figure(figsize = (10, 10))
    plt.plot(x, y, '-b', color='magenta', label='Hasil Pendekatan') 
    plt.plot(x, G, 'g--', color='blue', label='Hasil Analitik')
    plt.title(Judul)
    plt.xlabel('x')
    plt.ylabel('y = F(x)')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig('image\euler.png')
    print("hasil pendekatan yang diperoleh dari metode Euler adalah:", y)
    print("hasil analitik yang diperoleh adalah:", G)

#~~~METODE HEUN~~~
"""
Metode Heun adalah salah satu metode numerik yang dapat digunakan untuk menyelesaikan berbagai persoalan 
matematika yang mempunyai masalah nilai awal, yaitu masalah penyelesaian suatu persamaan diferensial dengan
syarat awal yang telah diketahui. Metode Heun juga merupakan salah satu peningkatan dari metode Euler. Metode ini melibatkan 
2 buah persamaan. Persamaan pertama disebut sebagai persamaan prediktor yang digunakan untuk memprediksi 
nilai integrasi awal dan persamaan kedua disebut sebagai persamaan korektor yang mengoreksi hasil integrasi awal.
"""
def Heun(h,x0,xn,y0,a,b,c,d):
    h = h
    x0 = x0
    xn = xn
    x = np.arange(x0, xn + h, h)
    y0 = y0
    a = a
    b = b
    c = c
    d = d
    G = a*(x**3) + b*(x**2) + c*x + d
    f = lambda x, y: a*(x**3) + b*(x**2) + c*x + d
    y = np.zeros(len(x))
    y[0] = y0

    for i in range(0, len(x) - 1):
        k1 = f(x[i], y[i])
        k2 = f((x[i]+h), (y[i]+(h*k1)))
        y[i+1] = y[i]+(0.5*h*(k1+k2))

    Galat = G-y
    print("galat yang diperoleh dari metode Heun adalah:", Galat)

    Judul = ("Grafik Pendekatan Persamaan Differensial Biasa Dengan Metode Heun")
    plt.figure(figsize = (10, 10))
    plt.plot(x, y, '-b', color='magenta', label='Hasil Pendekatan')
    plt.plot(x, G, 'g--', color='blue', label='Hasil Analitik')
    plt.title(Judul)
    plt.xlabel('x')
    plt.ylabel('y = F(x)')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig('image\heun.png')
    print("hasil pendekatan yang diperoleh dari metode Heun adalah:", y)
    print("hasil analitik yang diperoleh adalah:", G)

#~~~START~~~
print("Kode-kode modul: \n",
     "1. Modul Akar-Akar Persamaan \n",
     "2. Modul Sistem Persamaan Linier dan Matriks \n",
     "3. Modul Integrasi Numerik \n",
     "4. Modul Persamaan Diferensial Biasa \n")
setting = int(input("Masukkan kode penggunaan modul: "))
if (setting == 1):
    print("Kode-kode akar persamaan: \n",
          "1. Metode Setengah Interval \n",
          "2. Metode Interpolasi Linier \n",
          "3. Metode Newton-Rhapson \n",
          "4. Metode Secant \n",
          "5. Metode Iterasi \n")
    print("Persamaan; ax^3+bx^2+cx+d")
    setting = int(input("Masukkan kode penggunaan akar persamaan: "))
    if (setting == 1):
        X1 = float(input("Masukkan Nilai Pertama Setengah Interval: "))
        X2 = float(input("Masukkan Nilai Kedua Setengah Interval: "))
        a = float(input("Masukkan Nilai a Setengah Interval: "))
        b = float(input("Masukkan Nilai b Setengah Interval: "))
        c = float(input("Masukkan Nilai c Setengah Interval: "))
        d = float(input("Masukkan Nilai d Setengah Interval: "))
        X = Setengah_Interval(X1,X2,a,b,c,d)
        print(X)
    elif (setting == 2):
        X1 = float(input("Masukkan Nilai Pertama Interpolasi Linier: "))
        a = float(input("Masukkan Nilai a Interpolasi Linier: "))
        b = float(input("Masukkan Nilai b Interpolasi Linier: "))
        c = float(input("Masukkan Nilai c Interpolasi Linier: "))
        d = float(input("Masukkan Nilai d Interpolasi Linier: "))
        X = Interpolasi_Linier(X1,a,b,c,d)
        print(X)
    elif (setting == 3):
        X1 = float(input("Masukkan Nilai Pertama Newton-Rhapson: "))
        a = float(input("Masukkan Nilai a Newton-Rhapson: "))
        b = float(input("Masukkan Nilai b Newton-Rhapson: "))
        c = float(input("Masukkan Nilai c Newton-Rhapson: "))
        d = float(input("Masukkan Nilai d Newton-Rhapson: "))
        X = Newton_Rhapson(X1,a,b,c,d)
        print(X)
    elif (setting == 4):
        X1 = float(input("Masukkan Nilai Pertama Secant: "))
        a = float(input("Masukkan Nilai a Secant: "))
        b = float(input("Masukkan Nilai b Secant: "))
        c = float(input("Masukkan Nilai c Secant: "))
        d = float(input("Masukkan Nilai d Secant: "))
        X = Secant(X1,a,b,c,d)
        print(X)
    else:
        X1 = float(input("Masukkan Nilai Pertama Iterasi: "))
        a = float(input("Masukkan Nilai a Iterasi: "))
        b = float(input("Masukkan Nilai b Iterasi: "))
        c = float(input("Masukkan Nilai c Iterasi: "))
        d = float(input("Masukkan Nilai d Iterasi: "))
        X = Iterasi(X1,a,b,c,d)
        print(X)
elif (setting == 2):
    print("Kode-kode sistem persamaan linier dan matriks: \n",
          "1. Metode Gauss \n",
          "2. Metode Gauss-Jordan \n",
          "3. Metode Gauss-Siedel \n",
          "4. Metode Jacobi")
    setting = int(input("Masukkan kode penggunaan sistem persamaan linier dan matriks: "))
    if (setting == 1):
        print("Matriks: \n",
              "[[a1,a2,a3], [b1,b2,b3], [c1,c2,c3]] \n",
              "[D1,D2,D3]")
        a1 = int(input("Masukkan nilai a1:"))
        a2 = int(input("Masukkan nilai a2: "))
        a3 = int(input("Masukkan nilai a3: "))
        b1 = int(input("Masukkan nilai b1: "))
        b2 = int(input("Masukkan nilai b2: "))
        b3 = int(input("Masukkan nilai b3: "))
        c1 = int(input("Masukkan nilai c1: "))
        c2 = int(input("Masukkan nilai c2: "))
        c3 = int(input("Masukkan nilai c3: "))
        D1 = int(input("Masukkan nilai D1: "))
        D2 = int(input("Masukkan nilai D2: "))
        D3 = int(input("Masukkan nilai D3: "))
        A = np.array([[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]], dtype=int)
        B = np.array([D1, D2, D3])
        AB = np.hstack([A, B.reshape(-1, 1)])
        X = Gauss(A, B)
        print(X)
    elif (setting == 2):
        print("Matriks: \n",
              "[[a1,a2,a3,D1], \n",
              "[b1,b2,b3, D2], \n",
              "[c1,c2,c3, D3]]")
        n = int(input("Masukkan banyaknya jumlah variabel yang akan dicari: "))
        a1 = int(input("Masukkan nilai a1:"))
        a2 = int(input("Masukkan nilai a2: "))
        a3 = int(input("Masukkan nilai a3: "))
        b1 = int(input("Masukkan nilai b1: "))
        b2 = int(input("Masukkan nilai b2: "))
        b3 = int(input("Masukkan nilai b3: "))
        c1 = int(input("Masukkan nilai c1: "))
        c2 = int(input("Masukkan nilai c2: "))
        c3 = int(input("Masukkan nilai c3: "))
        D1 = int(input("Masukkan nilai D1: "))
        D2 = int(input("Masukkan nilai D2: "))
        D3 = int(input("Masukkan nilai D3: "))
        m = np.array([[a1, a2, a3, D1],
                      [b1, b2, b3, D2],
                      [c1, c2, c3, D3]],dtype=int)
        print('Matriks Persamaan')
        print(m)
        m = GaussJordan(m,n)
        print(f"""Hasil Pengelolaan menggunkan metode Gauss Jordan didapatkan matriks:
        {m}""")
    elif (setting == 3):
        print("Persamaan: \n",
              "x1=(Dx1-bx1-cx3)/ax1 \n",
              "x2=(Dx2-bx2-cx2)/bx2 \n",
              "x3=(Dx3-ax3-bx3)/cx3")
        a1 = float(input("Masukkan nilai a1: "))
        a2 = float(input("Masukkan nilai a2: "))
        a3 = float(input("Masukkan nilai a3: "))
        b1 = float(input("Masukkan nilai b1: "))
        b2 = float(input("Masukkan nilai b2: "))
        b3 = float(input("Masukkan nilai b3: "))
        c1 = float(input("Masukkan nilai c1: "))
        c2 = float(input("Masukkan nilai c2: "))
        c3 = float(input("Masukkan nilai c3: "))
        D1 = float(input("Masukkan nilai D1: "))
        D2 = float(input("Masukkan nilai D2: "))
        D3 = float(input("Masukkan nilai D3: "))
        r = int(input("Masukkan range iterasi: "))
        X = Gauss_Siedel(a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3,r)
        print(X)
    elif (setting == 4):
        print("Persamaan: \n",
              "x1=(Dx1-bx1-cx3)/ax1 \n",
              "x2=(Dx2-bx2-cx2)/bx2 \n",
              "x3=(Dx3-ax3-bx3)/cx3")
        a1 = float(input("Masukkan nilai a1 Jacobi: "))
        a2 = float(input("Masukkan nilai a2 Jacobi: "))
        a3 = float(input("Masukkan nilai a3 Jacobi: "))
        b1 = float(input("Masukkan nilai b1 Jacobi: "))
        b2 = float(input("Masukkan nilai b2 Jacobi: "))
        b3 = float(input("Masukkan nilai b3 Jacobi: "))
        c1 = float(input("Masukkan nilai c1 Jacobi: "))
        c2 = float(input("Masukkan nilai c2 Jacobi: "))
        c3 = float(input("Masukkan nilai c3 Jacobi: "))
        D1 = float(input("Masukkan nilai D1 Jacobi: "))
        D2 = float(input("Masukkan nilai D2 Jacobi: "))
        D3 = float(input("Masukkan nilai D3 Jacobi: "))
        r = int(input("Masukkan range iterasi: "))
        X = Jacobi(a1,a2,a3,b1,b2,b3,c1,c2,c3,D1,D2,D3,r)
        print(X)
elif (setting == 3):
    print("Kode-kode integrasi numerik: \n",
          "1. Metode Trapesium Satu Pias \n",
          "2. Metode Trapesium Banyak Pias \n",
          "3. Metode Simpson 1/3 \n",
          "4. Metode Simpson 3/8")
    print("Persamaan; ax^3+bx^2+c")
    setting = int(input("Masukkan kode penggunaan integrasi numerik: "))
    if (setting == 1):
        A = int(input("Masukkan Nilai Batas Bawah Integral: "))
        B = int(input("Masukkan Nilai Batas Atas Integral: "))
        a = int(input("Masukkan Nilai a: "))
        b = int(input("Masukkan Nilai b: "))
        c = int(input("Masukkan Nilai c: "))
        X = Trapesium_SatuPias(A,B,a,b,c)
        print(X)
    elif (setting == 2):
        A = int(input("Masukkan Nilai Batas Bawah Integral: "))
        B = int(input("Masukkan Nilai Batas Atas Integral: "))
        N = int(input("Masukkan Jumlah Pias: "))
        a = int(input("Masukkan Nilai a: "))
        b = int(input("Masukkan Nilai b: "))
        c = int(input("Masukkan Nilai c: "))
        X = Trapesium_BanyakPias(A,B,N,a,b,c)
        print(X)
    elif (setting == 3):
        A = int(input("Masukkan Nilai Batas Bawah Integral: "))
        B = int(input("Masukkan Nilai Batas Atas Integral: "))
        a = int(input("Masukkan Nilai a: "))
        b = int(input("Masukkan Nilai b: "))
        c = int(input("Masukkan Nilai c: "))
        X = Simpson_1per3(A,B,a,b,c)
        print(X)
    else:
        A = int(input("Masukkan Nilai Batas Bawah Integral: "))
        B = int(input("Masukkan Nilai Batas Atas Integral: "))
        a = int(input("Masukkan Nilai a: "))
        b = int(input("Masukkan Nilai b: "))
        c = int(input("Masukkan Nilai c: "))
        X = Simpson_3per8(A,B,a,b,c)
        print(X)
else:
    print("Kode-kode persamaan diferensial biasa: \n",
          "1. Metode Euler \n",
          "2. Metode Heun \n")
    print("Persamaan; ax^3+bx^2+cx^3+d")
    setting = int(input("Masukkan kode penggunaan persamaan diferensial biasa: "))
    if (setting == 1):
        h = float(input("Masukkan Nilai h Euler: "))
        x0 = float(input("Masukkan Nilai x0 Euler: "))
        xn = float(input("Masukkan Nilai xn Euler: "))
        y0 = float(input("Masukkan Nilai y0 Euler: "))
        a = float(input("Masukkan Nilai a Euler: "))
        b = float(input("Masukkan Nilai b Euler: "))
        c = float(input("Masukkan Nilai c Euler: "))
        d = float(input("Masukkan Nilai d Euler: "))
        X = Euler(h,x0,xn,y0,a,b,c,d)
        print(X)
    else:
        h = float(input("Masukkan Nilai h Heun: "))
        x0 = float(input("Masukkan Nilai x0 Heun: "))
        xn = float(input("Masukkan Nilai xn Heun: "))
        y0 = float(input("Masukkan Nilai y0 Heun: "))
        a = float(input("Masukkan Nilai a Heun: "))
        b = float(input("Masukkan Nilai b Heun: "))
        c = float(input("Masukkan Nilai c Heun: "))
        d = float(input("Masukkan Nilai d Heun: "))
        X = Heun(h,x0,xn,y0,a,b,c,d)
        print(X)


# 
