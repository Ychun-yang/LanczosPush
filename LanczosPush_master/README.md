# Theoretically and Practically Efficient Resistance Distance Computation on Large Graphs


## Compile:
```
sh runscript.sh
```

## Parameters:  
- -d \<datasets\> (e.g., dblp)
- -algo \<algorithm\> (e.g., lanczos, lanczos_push)
- -n \<number of query nodes\> 
- -e \<the Absolute Error parameter $\epsilon$\> 
- -k \<iteration number\>



## Remarks:
- The "data" directory is used for storing the dataset. 
- The "result" directory is used for storing the output results of different algorithms. 
- The "plot_result" directory is used for storing the average runtime/error output by different algorithms under different parameter settings. 
