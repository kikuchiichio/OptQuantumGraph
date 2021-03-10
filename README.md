# OptQuantumGraph
Computer programs used in "Optimization of quantum graphs using heuristic algorithms (Ichio Kikuchi, Akihito Kikuchi,March, 2021)"

These programs are made for following computations which are explained in the article. 

They require Lapack and Blas libraries.

(1) Protptype.c : Computation of the determinant to obtain the energy spectrum of the star-shaped graph.

(2) ppp.c : Optimization on nodes of the star-shaped graph.

(3) pppp.c : Optimization on edges of the star-shaped graph. 


These  three programs issue following lines through the computations.

### The number of the iteeration
 No. 135
### These lines present the index of the sampling particle, two-dimensional coorinates (in the left five columns),

### two-dimensional velocities, and the corresponding value of the objective functions (in the rightmost column).

0 0.942332 -0.443396 -0.000085 0.000657 0.691388  

1 0.942321 -0.443313 -0.000064 0.000496 0.691383

2 0.942347 -0.443511 -0.000058 0.000447 0.691396

3 0.942341 -0.443462 -0.000042 0.000329 0.691393

4 0.942326 -0.443321 -0.000042 0.000349 0.691384

5 0.942328 -0.443359 -0.000040 0.000314 0.691386

6 0.942302 -0.443189 -0.000075 0.000562 0.691375

7 0.942335 -0.443424 -0.000040 0.000325 0.691390

8 0.942322 -0.443327 -0.000049 0.000392 0.691384

9 0.942339 -0.443455 -0.000045 0.000351 0.691392

### This line presents the statics: 

### the averages and the variations in two-dimensional coorinates of the sampling particles (AVR and VAR); 

### the coordinates of the ever-best particle (BEST), 

### the corresponding vallue of the objective function (gbestval)

### and the determinant (det^2). 

STATISTICS: AVR=  0.942329  -0.443376  VAR= 0.000000  0.000000 BEST: 0.942302 -0.443189 gbestval=0.691375 det^2=0.00000




(4) ga.c : a discrete problem concerning the edge-connection.

You should run this program through 'gav.sh' with the input file ( simple_ga_input.txt).

The input-file assigns the range of the genes. You should not modify it.

The program issues the lines as follows. 

TARGET:0.3 BEST GENE: ( 1, 1, 1, 1, 1,) BEST VALUE:0.465

TARGET:0.4 BEST GENE: ( 1, 1, 1, 1, 1,) BEST VALUE:0.465

TARGET:0.5 BEST GENE: ( 1, 1, 1, 1, 1,) BEST VALUE:0.465

TARGET:0.6 BEST GENE: ( 0, 1, 1, 1, 1,) BEST VALUE:0.595

Each line show the TARGET VALUE, THE BEST GENE, and THE ACHIEVED VALUE OF THE BEST GENE.




March 10th, 2021,cWe are still ameliorating the contents. 
If you need them urgently or if you find some questionable places,
please send e-mails to the corresponding author (Akihito Kikuchi) of the article.
