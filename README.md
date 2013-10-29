Overview
--------

SEEDOFLIFE is a matlab toolbox for computing the seed-of-life nested shape descriptor.

This toolbox is released in conjunction with the ICCV'13 paper "[Nested Shape Descriptors](https://www.dropbox.com/s/g8yc76ffx8ia99d/iccv13_nsd_final.pdf)",
and can be used to reproduce results.

Installation
------------

You must run set_paths to initialize dependencies before running demos.

sh\> git clone https://github.com/jebyrne/seedoflife  
\>\> cd seedoflife  
\>\> set_paths  
\>\> run_iccv13  


Abstract
----------
In this paper, we propose a new family of binary local feature descriptors called nested shape descriptors.  These descriptors are constructed by pooling oriented gradients over a large geometric structure called the Hawaiian earring, which is constructed with a nested correlation structure that enables a new robust local distance function called the nesting distance.  This distance function is unique to the nested descriptor and provides robustness to outliers from order statistics.  In this paper, we define the nested shape descriptor family and introduce a specific family member called the seed-of-life descriptor.  We perform a trade study to determine optimal descriptor parameters for the task of image matching.  Finally, we evaluate performance compared to state-of-the-art local feature descriptors on the VGG-Affine image matching benchmark, showing significant performance gains.  Our descriptor is the first binary descriptor to outperform SIFT on this benchmark.    


Reference
----------
* J. Byrne and J. Shi, "Nested Shape Descriptors", International Conference on Computer Vision (ICCV'13), Sydney Australia, 2013

