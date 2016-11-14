# random_field
This is a python script designed to produce realistic-looking spatially correlated random field, such a s hydraulic conductivity, for use in 2-D or 3-D visualizations and/or numerical models. The model starts with a set of user specified seeds (locations with a known property), adds to the seed set sequentially by postulating new nearby points chosen from a Gaussian distribution, and then generates a numerical grid using scipy’s gridding routine once the seed population maximum is reached. The code is not mean to be geostatistically robust but rather an easy-to-understand demo that produces reasonable looking results.

The following tab-delimited input files are required (assisted with a tkinter-based GUI):

* domain.txt - grid settings
* params.txt - interpolation factors, seed generation
* seeds.txt - initial seed points (need at least one)

More background information can be found here: https://numericalenvironmental.wordpress.com/2016/07/11/random-field-generation/

An example application can be found here: https://numericalenvironmental.wordpress.com/2016/10/03/napl-migration-through-a-correlated-random-field/

Email me with questions at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
