# Range Queries Code (Under Construction)

This repository includes experiments from the  paper 
> [Answering Range Queries Under Local Differential Privacy](https://arxiv.org/abs/1812.10942)

# The experiments involve 3 histogram aggregation primitives. 
+ OUE ([Optimal Unary Encoding](https://www.usenix.org/system/files/conference/usenixsecurity17/sec17-wang-tianhao.pdf))
+ OLH ([Optimal Local Hashing](https://www.usenix.org/system/files/conference/usenixsecurity17/sec17-wang-tianhao.pdf)) 
+ HRR ([Hadamard Randomized Response](http://dimacs.rutgers.edu/~graham/pubs/papers/sigmod18.pdf)) 


## Dependencies
+ [CMake](https://cmake.org/download/) (>= 2.8.7) 
+ C++11 support required 

## Repository Structure 
+ Headers (.h header files.)
+ sources (.cpp implementation files.)
+ Main (Driver files implementing the experiments from the paper.)
+ treefiles (This directory contains cache for the Haar transform based method (**HaarHRR**).
 Without these files, the HaarHRR method may take a long time when the D is large. Please make sure that the treefiles directory is downloaded.)
+ Misc (Misc. cpp and headers implementing HRRk (sampling and perturbing k coeffcients instead of one coefficient.)) This directory can be ignored.

 The code has been tested on a Linux machine. I haven't tested thoroughly, but I hope it can be easily ported on other planforms.

 
## Main files description.
 
 * mseranges.cpp -- The function **mseranges()** computes the mean squared error for all interval ranges lengths for various fan outs. 
 * avgerr.cpp -- The function **avg_error_eps()** compute the over all mean squared error for a various fan outs and epsilon values. The second function 
 **avg_error_dbn()** also computes the over all mean squared error for a various fan outs but for various distribution shapes parameterized by the value **prob**. 

## Compilation steps.

 Step 1, clone the repository.
```bash
    $ git clone https://gitlab.com/Tejasvk/RangeQueries
```

 Step 2
 
 Build the project. Execute

 ```bash
 $ cd RangeQueries/
 $ cmake .
 $ make
 ```
 
This step builds the static library **librangequeries.a**.

 Step 3
     Now compile the code:
 ```bash    
    $ cd main
    $ g++ -std=c++11 -O3 -o avgerror.cpp avgerror.out librangequeries.a
    $ ./avgerror.out
 ```

## To Do
* Some more refactoring.
* Add more plotting scripts.
* Add more documentation.
* Add tests.

Please include a citation to the paper if you use this code.

>@inproceedings{CKS19,
  title = {Answering Range Queries Under Local Differential Privacy},
  author = {Graham Cormode and Tejas Kulkarni and Divesh Srivastava},  
  url  = {https://arxiv.org/abs/1812.10942},
  year = {2019},
}

Contact: [Tejas Kulkarni](https://warwick.ac.uk/fac/sci/dcs/people/research/u1554597/) at abc@gmail.com, where a=tejas, b=vijay, c=kulkarni 

