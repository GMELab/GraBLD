GraBLD

Gradient boosted and LD adjusted (GraBLD) is an R based software package that applied to polygenic traits prediction using gene scores. 

###Quick Start### 

The current release is: version 1.0. See "release" tab. 

To install our R package “GraBLD”, you can either run in R directly:

    install.packages("devtools") # if you have not installed already
    devtools::install_github("GMELab/GraBLD")

OR

After downloading GraBLD_0.1.0.tar.gz at "release" tab, run:

    R CMD INSTALL GraBLD_0.1.0.tar.gz

To load the library in R, simply run in R:

    library(“GraBLD”)

###Tips for using GraBLD on large datasets### 

The package contains two main functions, one is to calcualte the LD adjustments and the other to calculate the gradient boosted polygenic score weights. The two are combined in the end to give the actual weights of the polygenic score. 

##LD adjustments## 

For large datasets, it is recommended to run from the command line with 

    for((i = 1; i <= chr; i++))
    do
    Rscript PerformLDadj.R size data_name ${i} &
    done

where the R script "PerformLDadj.R" might look something like this, while additional options can be added to the argument list: 
    
    #!/bin/sh
    rm(list = ls())
    library('GraBLD')
    args = (commandArgs(TRUE))
    size = eval(parse(text=args[1]))
    source_data = args[2]
    chr = eval(parse(text=args[3]))
    geno_data = load_geno(source_data = source_data, PLINK = TRUE, chr = chr)
    geno_norm = full_normal_geno(geno_data)
    LD_OUT <- LDadj(geno_raw = geno_norm, chr = chr, size = size, write = TRUE)
  
  
##Gradient boosted weights## 

For large datasets, it is recommended to run from the command line with 
    
    geno_data="chr1.raw" # or the actual file name of the genotype data 
    trait_name="BMI" # name of the trait
    annotations_file="annotation.txt" # or the actual file name of the annotation data
    validation=5 # or value of your choosing
    interaction_depth=5 # or value of your choosing
    shrinkage_parameter=0.01 # or value of your choosing
    bag_fraction=0.5 # or value of your choosing
    maximum_tree=2000 # or value of your choosing
    for (( i = 1; i <= $validation; i++))
    do
    Rscript calculate_gbm.R $geno_data
       $trait_name $annotations_file ${i}
       $validation $interaction_depth
       $shrinkage_parameter $bag_fraction
       $maximum_tree &
    done
    
where the R script calculate_gbm.R might look something like this, while additional options can be added to the argument list: 

    #!/bin/sh
    rm(list = ls())
    library('GraBLD')
    args = (commandArgs(TRUE))
    geno_data = args[1]
    trait_name = args[2]
    annotations_file = args[3]
    steps = eval(parse(text=args[4]))
    validation = eval(parse(text=args[5]))
    p1 = eval(parse(text=args[6]))
    p2 = eval(parse(text=args[7]))
    p3 = eval(parse(text=args[8]))
    p4 = eval(parse(text=args[9]))
    betas = load_beta(trait_name)
    annotation = load_database(annotations_file, pos = 2:3) # taking the 2nd and 3rd columns of "annotations_file"
    geno <- load_geno(geno_data)
    GraB(betas = betas, annotations = annotation,
      trait_name = trait_name, steps = steps, validation = validation,
      interval = 200, sig = 1e-05, interact_depth = p1, shrink = p2,
      bag_frac = p3, max_tree = p4, WRITE = TRUE)


##Combining the two## 

The two steps above produce separate files, one for the LD adjustments and the other the boosted weights, and they can be combined by the function "GraBLD.score". 

    LD_val <- read.table("OUTPUTS_FROM_STEP1.txt") ## if each chromosome is computed separately, you will need to append them into a single file 
    gbm_val <- read.table("OUTPUTS_FROM_STEP2.txt")
    gs <- GraBLD.score(geno_raw = YOUR_GENO_DATA, LDadjVal = LD_val, gbmVal = gbm_val, Pheno = NULL)

Notice that if "Pheno" is supplied, both the polygenic gene score as well as the prediction R-squared (adjusted) are returned, otherwise only the polygenic gene score is returned.
