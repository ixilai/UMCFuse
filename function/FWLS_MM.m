%  FGS  Fast global image smoother.
% 
%  F = FGS(img, sigma, lambda, joint_image, num_iterations, attenuation)
%
%  Parameters:
%    img             Input image to be filtered [0,255].
%    sigma           Filter range standard deviation.
%    lambda          Filter lambda.
%    joint_image     (optional) Guidance image for joint filtering [0,255].
%    num_iterations  (optional) Number of iterations to perform (default: 3).
%    attenuation     (optional) attenuation factor for iteration (default: 4).
%
%  This is the reference implementation of the fast global image smoother
%  described in the paper:
% 
%    Fast Global Image Smoothing based on Weighted Least Squares
%    D. Min, S. Choi, J. Lu, B. Ham, K. Sohn, and M. N. Do,
%    IEEE Trans. Image Processing, vol. no. pp., 2014.
%
%  Please refer to the publication above if you use this software. For an
%  up-to-date version go to:
%  
%  https://sites.google.com/site/globalsmoothing/
%
%  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY EXPRESSED OR IMPLIED WARRANTIES
%  OF ANY KIND, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
%  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%  OUT OF OR IN CONNECTION WITH THIS SOFTWARE OR THE USE OR OTHER DEALINGS IN
%  THIS SOFTWARE.
%
%  Version 1.0 - December 2014.

function F = FWLS_MM(img, img_g, wtbx, wtby, lambda, num_iterations, attenuation)

    I = double(img);

    if ~exist('num_iterations', 'var')
        num_iterations = 3;
    end
    
    if ~exist('attenuation', 'var')
        attenuation = 4;
    end
    %wtbx= bilateralFilter(wtbx, wtbx, min(wtbx(:)), max(wtbx(:)), 3, 0.5);
    
    F = mexFWLS_MM(I, img_g, wtbx, wtby, lambda, num_iterations, attenuation);
    
    %% The intuitive mex version
%     F = mexFGS_simple(I, J, sigma, lambda, num_iterations, attenuation);
    
    %% Return the result
    F = cast(F, class(img));
    
end