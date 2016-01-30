# Network-ILP

Network optimization with integer linear programming.

If you are using this code please cite:

```
@article{Rempfler2015,
	author = "Markus Rempfler and Matthias Schneider and Giovanna D. Ielacqua and Xianghui Xiao and Stuart R. Stock and Jan Klohs and Gábor Székely and Bjoern Andres and Bjoern H. Menze",
	title = "Reconstructing cerebrovascular networks under local physiological constraints by integer programming",
	journal = "Medical Image Analysis ",
	volume = "25",
	number = "1",
	pages = "86 - 94",
	year = "2015",
	issn = "1361-8415",
}
```


### Requirements

CMake

Boost

IBM ILOG CPLEX (non-free, **not included**)

### Installation 

Starting in the top level directory:

```
mkdir build
cd build
cmake ../src
make install
```

### License

Copyright (c) 2015, Markus Rempfler and Bjoern Andres.

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* The name of the author must not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
