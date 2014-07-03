##g_resid_distrib
***

###About

This program calculates the positions of center of masses of the selected
groups along X, Y and Z axes. This tool is useful to calculate the position
of selected residues along the channel-axis of the protein channel. Further,
output files could be use to obtain the distribution of residues' position
during the MD simulations.

***

###Requirements
To compile and install, GROMACS libraries <code> libgmx, libmd, libgmxana </code> are required.
***

###Download
<pre><code>git clone https://github.com/rjdkmr/g_resid_distrib
</code></pre>
***

###Installation
<pre><code>cd g_resid_distrib
mkdir build
cd build
cmake ..  -DGMX_PATH=/opt/gromacs -DCMAKE_INSTALL_PREFIX=/opt/g_resid_distrib
make
make install
</code></pre>

Directory <code>/opt/gromacs</code> should contains <code>include</code> and <code> lib </code> directories. If these directories are in seprate locations, use followings:
<pre><code>cmake ..  -DGMX_LIB=/path/to/lib -GMX_INCLUDE=/path/to/include -DCMAKE_INSTALL_PREFIX=/opt/g_resid_distrib
</code></pre>

If fftw library <code> libfftw3f.so or libfftw3f.a </code> are not present in standard locations:
<pre><code>-DFFTW_LIB=/path/to/fftw3/lib</code></pre>
***

###Usage
<pre><code>g_resid_distrib -h
</code></pre>
***
