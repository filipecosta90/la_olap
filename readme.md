## Welcome to LaOLAP - a Linear Algebra database engine

### What is this?

This is a Multi-Core Linear Algebra based database engine written to use vector based processing and in an near future high bandwidth of modern GPUs and Co-processors.

### Requirements
* bison
* flex
* icc 
* Intel Math Kernel Library 

### How to build?

```
git clone https://github.com/filipecosta90/la_olap.git
cd la_olap/src
make
```

### Will Feature :

  *  Vector-based processing  
	
  * Column-based storage  
    Minimizes disk I/O by only accessing the relevant data.

