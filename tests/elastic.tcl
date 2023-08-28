# # t     p         ph       a        v       d       u
#  0.0  0.0000             0.0000  0.0000  0.0000  0.0000 
#  0.1  5.0000    5.0000  17.4666  0.8733  0.0437  0.0328 
#  0.2  8.6603   26.6355  23.1801  2.9057  0.2326  0.2332 
#  0.3 10.0000   70.0837  12.3719  4.6833  0.6121  0.6487 
#  0.4  8.6603  123.9535 −11.5175  4.7260  1.0825  1.1605 
#  0.5  5.0000  163.8469 −38.1611  2.2421  1.4309  1.5241 
#  0.6  0.0000  162.9448 −54.6722 −2.3996  1.4230  1.4814 
#  0.7  0.0000  110.1710 −33.6997 −6.8182  0.9622  0.9245 
#  0.8  0.0000   21.8458  −2.1211 −8.6092  0.1908  0.0593 
#  0.9  0.0000  −69.1988  28.4423 −7.2932 −0.6043 −0.7751 
#  1.0  0.0000 −131.0066  47.3701 −3.5026 −1.1441 −1.2718

# set period 1.0
# set M [expr 1.0];
# set K [= 4.0*pi*pi*$M/($period*$period)];
# # set zeta   0.001
# # set C [expr 2.0*sqrt(K/M)*zeta];

set K  10.0;
set M  0.2533
set C  0.1592
set uy 0.75
set options "-beta [= 1/6] -gamma 0.5" ; # linear acceleration
# set options "-beta [= 1/4] -gamma 0.5" ; # constant acceleration

set dt 0.1
set n 100
# set series  [lmap i [linspace 0 10 $n] {expr sin($i)}]
set series {
 0.0000
 5.0000
 8.6603
10.0000
 8.6603
 5.0000
 0.0000
 0.0000
 0.0000
 0.0000
 0.0000
}

model basic 1 1

uniaxialMaterial Elastic 1 $K
# uniaxialMaterial ElasticPP 1 $K $uy


invoke UniaxialMaterial 1 {
  puts [
    string map {" " "\n"} [integrate $dt $series -mass $M -damp $C {*}$options]
  ]
}

