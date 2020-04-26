## Agile Combinatorial Factor Decomposition (AgileFD)

This code accompanies the following publications:

Bai, J., Xue, Y., Bjorck, J., Le Bras, R., Rappazzo, B., Bernstein, R., Suram, S.K., Van Dover, R.B., Gregoire, J.M. and Gomes, C.P., 2018. ***Phase mapper: Accelerating materials discovery with AI.*** *AI Magazine*, *39*(1), pp.15-26.

Xue, Y., Bai, J., Le Bras, R., Rappazzo, B., Bernstein, R., Bjorck, J., Longpre, L., Suram, S., van Dover, B., Gregoire, J., Gomes, C (2017). ***Phase-Mapper: An AI Platform to Accelerate High Throughput Materials Discovery.*** In Twenty-ninth AAAI Conference on Artificial Intelligence (AAAI'17).

Please cite these or related works as appropriate in follow-up publications.

Bibtex:

```
@article{bai2018phase,
  title={Phase mapper: Accelerating materials discovery with AI},
  author={Bai, Junwen and Xue, Yexiang and Bjorck, Johan and Le Bras, Ronan and Rappazzo, Brendan and Bernstein, Richard and Suram, Santosh K and Van Dover, Robert Bruce and Gregoire, John M and Gomes, Carla P},
  journal={AI Magazine},
  volume={39},
  number={1},
  pages={15--26},
  year={2018}
}
@inproceedings{xue2017phase,
  title={Phase-Mapper: an AI platform to accelerate high throughput materials discovery},
  author={Xue, Yexiang and Bai, Junwen and Le Bras, Ronan and Rappazzo, Brendan and Bernstein, Richard and Bjorck, Johan and Longpre, Liane and Suram, Santosh K and van Dover, Robert B and Gregoire, John and others},
  booktitle={Twenty-Ninth IAAI Conference},
  year={2017}
}
```



### Introduction and Prerequisites

AgileFD is an Non-negative Matrix Factorization (NMF) based decomposition method. It is simple, agile and easy to implement. AgileFD is the prototype of many follow-up variants like IAFD, PIAFD, etc.

AgileFD has two external library dependencies which must be installed:

1) Armadillo, see: http://arma.sourceforge.net
2) ILOG/CPLEX Optimization Studio: https://www.ibm.com/developerworks/downloads/ws/ilogcplex/

AgileFD also uses TCLAP (http://tclap.sourceforge.net/) to process command-line arguments and Tino Kluge's spline implementation (http://kluge.in-chemnitz.de/opensource/spline/). These are header-only includes and are provided in the source distribution.

### Compilation

An example Makefile for GNU make and g++ is included, but will need to be modified to reflect the installation location of CPLEX, and the BLAS library (like openblas) used with Armadillo. Additional changes will be needed for alternate compilers. This configuration has been used on a few different Linux distributions, but additional changes could be required for other systems.

After installing all the prerequisites, one should be able to compile the code with the following command:

```bash
make
```

If the compilation fails, try to `make clean` first.

### Usage

```bash
./agilefd [-h] [OPTIONS] --inst INSTANCE_FILENAME --m M --k K --sol SOLUTION_FILENAME
```

Options:

```
   --length <int>
     The customized length of the array of qvalues (could be different from
     the length of the given qvalues in the instance file). This value 
	 overrides --stepsize if both are specified.
   --qmax <double>
     The customized maximum q value (could be different from the maximum
     qvalue in the instance file)
   --qmin <double>
     The customized minimum q value (could be different from the minimum
     qvalue in the instance file)
   --snapshot <int>
     Every ? iterations, output one snapshot of the lastest solution
   --gibbs
     whether or not to enforce Gibbs phase rule
   --mipgap <double>
     mipgap for MIP done via CPLEX, the default is 0.1
   --sparsity <double>
     The overall sparsity coefficient
   --stepsize <double>
     Customized stepsize for resampling. A default stepsize of 0 will result
     in use of the standard stepsize. If both are specified, --length will take
	 precedence.
   --sampleInit <string>
     The filename containing initialzation from single-phase sample points
   --valueInit <string>
     The filename containing initialization values for phases and phase freezing
   --seed <int>
     Random seed for the random number generator. The default value of -1
     means the random seed will be time(0)).
   --c <double>
     Related to termination criterion: In one iteration, if
     (old_cost-new_cost)<c*old_cost, then the loop terminates
   --time <int>
     The maximum time(seconds) allowed to train the model
   --m <int>
     (required)  The number of possible different shifts
   --k <int>
     (required)  The number of phases
   --sol <string>
     (required)  The output filename of the solution
   --inst <string>
     (required)  Input instance file
   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.
   --version
     Displays version information and exits.
   -h,  --help
     Displays usage information and exits.
```



### Example Command

```bash
./agilefd --inst input/my_inst.txt --m 10 --k 8 --sol output/my_output.txt --time 1000 --c 1e-5 --sparsity 1.0 --gibbs --mipgap 0.1 --valueInit input/my_valueInit.txt --snapshot 300 --qmin 11.0 --qmax 59.0 --length 2000
```

The formats of the input instance file and the output solution file can be found in the following paper:

Le Bras, R., Bernstein, R., Gregoire, J. M., Suram, S. K., Gomes, C. P., Selman, B., & van Dover, R. B. (2014). A Computational Challenge Problem in Materials Discovery: Synthetic Problem Generator and Real-World Datasets. In Twenty-Eighth International Conference on Artificial Intelligence (AAAI 14).

[NOTE additions to solution file specific to AgileFD - H,HS,L]

### Appendix

#### Value Initialization File Format

A value initialization file specifies the basis patterns to use for initialization, as well as related configuration options. The format is as follows:

    // Comments can be included as complete lines beginning with two slashes
    
    // Q values corresponding to the basis vectors, in the same units as in the instance file (e.g. nm^-1 or A^-1). 
    // Basis vectors are resampled, so these values do not need to exactly match the ones in the instance file.
    
        Q=1,1.1,1.2,... 
    
    // Basis patterns are only shifted to the right (positive shift) in AgileFD, 
    // so initial or frozen basis patterns should be specified as far
    // to the left as expected. The V parameter is a multiplicative shift which 
    // is applied to the Q vector and affects all basis patterns in the same way.
    // The following effectively shifts all basis patterns 1% to the left.
    
        V=0.99
    
    // B1...BK specify the basis patterns you wish to initialize. The indices must be 
    // sequential, but fewer than K can be specified if desired.
    
        B1=0.3123,0.545234,......
        B2=0.4324,0.454243,......
        B3=0.42345,1.42344,......
        ...
    
    // Whether to freeze each basis pattern: 0=no, 1=yes
    
        F1=1
        F2=0
        F3=1
        ...
    
    // Lists sample indices in which each phase is allowed to appear. For example, 
    // phase 1 could appear at sample points 1,2,54,65,76......, etc. If S is not
    // specified for an initialized basis, it can appear at any sample by default.
    
        S1=1,2,54,65,76,......
        S2=1,2,4,7,111,......
        S3=67,89,111,......
    ...

#### Sample Initialization File Format

Sample initialization is an alternative to value initialization, where the desired basis 
patterns are taken from samples in the instance file. To initialize phase 1 using sample 
point 123 and phase 2 from sample point 56:

    B1=123
    B2=56

Q should not be specified because the patterns come from the instance file. However, 
V, F, and S are set the same as in a value initialization file.

### Contact me

If you have any questions or suggestions, please feel free to contact jb2467@cornell.edu or rab38@cornell.edu.
