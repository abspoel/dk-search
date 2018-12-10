This code can be used to search for primes p that have large values of d_k(p) and d_k^\*(p), as defined in our paper:

> Mark Abspoel and Niek J. Bouman and Berry Schoenmakers and Niels de Vreede. Fast Secure Comparison for Medium-Sized Integers and Its Application in Binarized Neural Networks. Cryptographers Track RSA Conference, 2019. To appear.

To use, install Rust using [rustup](https://rustup.rs/). Clone the repository and use `cargo run --release -- [args]` to execute.

Two helper scripts `run.py` and `table.py` are provided to start multiple processes to speed up the computation, and aggregate the resulting log files, respectively.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

