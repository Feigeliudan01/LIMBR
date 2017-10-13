# MSLearn
Applications of machine learning to proteomics data processing

# 
Because of the way docker build executes layers, it is important to feed it bite size chunks, unfortunately there is no good way to iterate this process of layering in docker build so it was necessary to hard code each operation as a separate RUN command.  The more straightforward options of shell based for loops or make files being executed by a run command results in hard to diagnose hangs during the docker build process, likely do to the large amount of data being moved in each layer.  This structure has the additional benefit of taking advantage of dockers native caching to speed the build process.