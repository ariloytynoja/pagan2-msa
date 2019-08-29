
[![PAGAN](http://wasabiapp.org/download/theme/icons/pagan.png)](http://wasabiapp.org/software/pagan/) 
# PAGAN2
PAGAN2 is a general-purpose method for the alignment of DNA, codon and amino-acid sequences as graphs. It aligns sequences either with pileup or, when related by a tree, using phylogeny-aware progressive alignment algorithm. In both cases it uses graphs to describe the uncertainty in the presence of characters at certain sequence positions. PAGAN2 is largely compatible with PAGAN but implements new algorithms for alignment anchoring and memory handling. PAGAN2 can align sequences of several hundreds of kilobases in length.

PAGAN2 uses the NCBI TOOLKIT and Boost libraries. The NCBI TOOLKIT does not compile with GCC compiler newer than 4.8 and also requires a rather old version of libc. The code compiles well on Ubuntu 14.04 and instructions are provided for building it with that using Docker.

If you **use Linux** and just want to **use PAGAN2** (i.e. not see the source code and compile it from scratch), please go to the PAGAN homepage at http://wasabiapp.org/software/pagan and download the latest PAGAN2 package there. This package includes all helper applications needed e.g. for the guidetree inference. If you **do not use Linux**, you can try to use PAGAN2 inside the precompiled Docker container available at Docker hub. Please read the instructions for that below.

### Instructions
Both options naturally require that one has Docker installed...

#### 1. Using the precompiled Docker image from Docker hub

**Get the image and rename it as pagan2**
```
docker pull ariloytynoja/pagan2
docker tag ariloytynoja/pagan2 pagan2
```

**Get some data**
```
wget https://raw.githubusercontent.com/ariloytynoja/pagan-msa/master/examples/454_pileup/454_reads.fas 
```

**Run PAGAN2 inside a container**
```
docker run --rm -v `pwd`:/data \
  pagan2 --pileup --homopolymer \
  -q 454_reads.fas -o 454_aligned
```
Here Docker the options tell to remove the container after finishing (``--rm``) and map the file systems (``-v `pwd`:/data``). As the options are always the same, we can write them in a script file and call that instead.

**Create a script file to simplify the command**
```
cat > pagan2.sh << EOF 
#!/bin/bash
docker run --rm -v `pwd`:/data pagan2 "\$@"
EOF

chmod +x pagan2.sh 

./pagan2.sh --pileup --homopolymer \
  -q 454_reads.fas -o 454_aligned
```

#### 2. Compiling PAGAN2 from scratch using Docker

**Get the dockerfiles**
```
 git clone https://github.com/ariloytynoja/pagan2-docker-build.git
 cd pagan2-docker-build/
 ```
 
**Build the libs and the app**
```
docker build \
    -t pagan2libs \
    -f Dockerfile.pagan2libs \
    .

docker build \
    -t pagan2app \
    -f Dockerfile.pagan2app \
    .
 ```
 
**Option 1: Get the statically linked binary out**
```
docker create \
    --name pagan2app \
    pagan2app

docker cp pagan2app:/pagan2/src/pagan2 ./pagan2
 ```
 The binary ```pagan2``` contains all dependencies and can be used on all modern Linux systems (Ubuntu 14.04 and later).
 
 
 **Option 2: Build the helper apps and get them out**
```
docker build \
    -t pagan2app \
    -f Dockerfile.pagan2apps \
   .
   
docker create \
    --name pagan2apps \
    pagan2apps

docker cp pagan2apps:/pagan2/pagan2apps.tgz ./pagan2apps.tgz

tar xzf pagan2apps.tgz

./pagan2/bin/pagan2 --pileup --homopolymer \
  -q 454_reads.fas -o 454_aligned
```

 **Option 3: Create a small PAGAN2 image and run it as above**
 ```
 docker build \
    -t pagan2 \
    -f Dockerfile.pagan2 \
    .
    
docker run --rm -v `pwd`:/data \
  pagan2 --pileup --homopolymer \
  -q 454_reads.fas -o 454_aligned
```
The Docker container should run on non-Linux systems.
