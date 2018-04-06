# Rectangular v1.6
Create orthogonal grids on the local Earth curve

## User interface
![alt text](https://github.com/larryyin/rectangular/blob/master/demo/demo.JPG?raw=true "User interface")

## Demo run
"GET /_draft?latO=40.74589374151832&lngO=-74.03184610270273&latA=40.73256186902815&lngA=-74.05031745520023&latB=40.75922266814361&lngB=-74.01336734381107&latC=40.728756956916236&lngC=-74.0202337988892&cs=10&abcPath= HTTP/1.1" 200 -  
DEM path: dem/example.tif  
Buildings path: buildings/hob.kml  
NLCD path: nlcd/nynj.tif  
Output path: out/  
Iteration 1: 0 x 0  
Iteration 2: 10.014300507948885 x 10.04115576345354  
Iteration 3: 10.007140026562181 x 10.041155763456038  
Iteration 4: 10.003564906598484 x 10.041155762298603  
Iteration 5: 10.00356492388977 x 10.020493245752231  
Iteration 6: 10.00356490603221 x 10.01012030057025  
Iteration 7: 10.003564910987649 x 10.005039382798607  
H1 DescribeResult(nobs=89614, minmax=(10.003552143154204, 10.003576663897945), mean=10.00356484120473, variance=3.0255225672309684e-11, skewness=-0.05156601662124588, kurtosis=-0.9292146387012457)  
H2 DescribeResult(nobs=89614, minmax=(10.002481115507411, 10.002536605873729), mean=10.002509299413159, variance=2.377449689285769e-10, skewness=-0.00838633393946396, kurtosis=-1.1952970306170112)  
ANG DescribeResult(nobs=89614, minmax=(-9.66668085665114, -9.642457893510567), mean=-9.65456529455321, variance=3.453701333226129e-05, skewness=-0.00031071674027270874, kurtosis=-1.0910253952128157)  
Creating output file that is 1908P x 1780L.  
Processing input file dem/example.tif.  
Using internal nodata values (e.g. -3.40282e+38) for image dem/example.tif.  
Copying nodata values from source dem/example.tif to destination tmp/dem_clipped.tif.  
0...10...20...30...40...50...60...70...80...90...100 - done.  
Extracting depths from DEM...  
DEM cellsize: 0.000020 x -0.000020  
Importing buildings...  
Imported from buildings/hob.kml: 438 buildings.  
Extracting NLCD classes from raster...  
NLCD raster cellsize: 0.000374 x -0.000371  
Write to output...  
Job completed successfully.

## Output
[Demo output](https://github.com/larryyin/rectangular/tree/master/demo/output "demo output")  

Stats
Nodes: 259 x 346 = 89614  
Extent: -74.050317, 40.728670, -74.013296, 40.763164  
Area: 8966844.08 m^2  
H1: mean(10.00), median(10.00), min(10.00), max(10.00)  
H2: mean(10.00), median(10.00), min(10.00), max(10.00)  
ANG: mean(-9.65), median(-9.65), min(-9.67), max(-9.64)  
Depth: mean(-10.070), median(-2.011), min(-74.743), max(19.518)  
Min time step along I: 0.376 s  
Min time step along J: 0.377 s  

![alt text](https://github.com/larryyin/rectangular/blob/master/demo/demo_output.JPG "Output demo")

# conda environment
conda env create -f environment.yml

# MIT License

Copyright (c) 2018 larryyin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
