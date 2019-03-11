
[![PAGAN](http://wasabiapp.org/download/theme/icons/pagan.png)](http://wasabiapp.org/software/pagan/) 
# PAGAN2
PAGAN2 is a general-purpose method for the alignment of DNA and amino-acid sequences as graphs. It aligns sequences either with pileup or, when related by a tree, using phylogeny-aware progressive alignment algorithm. In both cases it uses graphs to describe the uncertainty in the presence of characters at certain sequence positions. PAGAN2 is largely compatible with PAGAN but implements new algorithms for alignment anchoring and memory handling. PAGAN2 can align sequences of several hundreds of kilobases in length.

PAGAN2 uses the NCBI TOOLKIT and Boost libraries. The NCBI TOOLKIT does not compile with GCC compiler newer than 4.8 and also requires a rather old version of libc. The code compiles well on Ubuntu 14.04 and instructions are provided for building it with that using Docker.

If you just want to **use PAGAN2** (i.e. not see the source code and compile it from scratch), please download the statically linked binary from the "bin/"" directory. The binary should work on all modern Linux versions. If you are not using Linux, you can try to use PAGAN2 inside the precompiled Docker container avaialbel at Docker hub.

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
mkdir /home/$USER/data
wget https://raw.githubusercontent.com/ariloytynoja/pagan-msa/master/examples/454_pileup/454_reads.fas -P /home/$USER/data
```

As we move data between our own file system and the Docker container, we need to map a directory in the first to a directory in the second. For simplicity we create a directory called "data" in our home directory and map that to a directory called "/data" in the container. 

**Run PAGAN2 inside a container**
```
docker run --rm -v /home/$USER/data:/data \
  pagan2 --ncbi --pileup --homopolymer \
  -q /data/454_reads.fas -o /data/454_aligned
```
Here the options tell to remove the container after finishing (``--rm``) and map the file systems (``-v /home/$USER/data:/data``). Note that the input and output filepaths have to start with "/data/"!

As the options are always the same, we can write them in a script file and call that instead.

**Creating a script file to simplify the command**
```
cat > pagan2_docker << EOF 
#!/bin/bash
docker run --rm -v /home/$USER/data:/data pagan2 "\$@"
EOF

chmod +x pagan2_docker 

./pagan2_docker --ncbi --pileup --homopolymer \
  -q /data/454_reads.fas -o /data/454_aligned
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
 
**Get the statically linked binary out**
```
docker create \
    --name pagan2app \
    pagan2app

docker cp pagan2app:/pagan2/src/pagan2 ./pagan2
 ```
 The binary ```pagan2``` contains all dependencies and can be used on all modern Linux systems (Ubuntu 14.04 and later).
 
 **Create PAGAN2 app image and run it**
 ```
 docker build \
    -t pagan2 \
    -f Dockerfile.pagan2 \
    .
    
docker run --rm -v /home/$USER/data:/data \
  pagan2 --ncbi --pileup --homopolymer \
  -q /data/454_reads.fas -o /data/454_aligned
```
The Docker container should run on non-Linux systems.
