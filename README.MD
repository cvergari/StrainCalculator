# StrainCalculator

StrainCalculator calculates strains of a flat 2D regular grid 
undergoing a series of deformation steps, resulting in displacements of
the original grid. A typical application is the post-processing of
Digital Image Correlation (DIC) methods. In brief, a grid is somehow
applied to a sample, the sample is then loaded in a given number of
loading steps, and the position of each grid element is retrieved at
each step. An example of grid deformation:

![Series of grid deformations](docs/_images/pointcloud.gif)
*An example of a planar grid undergoing a series of deformations. Source: [µDIC](https://mudic.readthedocs.io/en/latest/index.html)* 



## Context

Several toolboxes are available online for the processing of DIC data and calculation of strains. However, the strain calculation is often well hidden or not well documentent. The aim of this toolbox is to propose different methods to calculate strains, with clear documentation and examples.
This program is aimed at engineers, engineering students, but also researchers who use DIC and want to better control what is being calculated.



## Usage 

StrainCalculator() can use the displacements of the grid at each step
to compute a strain map.

```matlab

    obj = StrainCalculator();
    obj.generateFakeData('HorizontalCompression');
    obj.calculateStrains;
      
    for k = 1 : obj.nSteps
        cla
        obj.plotStrain(k);
        pause(0.3)
	end

```
	
Strains should be provided as an X and Y array of cells, each
containing the X and Y positions of the grid elements at each time
step.
For more detailed instructions look at the [documentation](docs/index.html) or directly at the [Getting Started](docs/GettingStarted.html) section.


# Version [0.1]


## License

MIT License

Copyright (c) [2022] [Claudio Vergari]

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

[MIT](https://choosealicense.com/licenses/mit/)


