# Cryptanalysis of the RSA variant based on cubic Pell equation

## Introduction

This is a Python implementation of lattice-based attack proposed in **Cryptanalysis of the RSA variant based on cubic Pell equation**[^cPRSA]. Some underlying functions are based on [Joachim Vandersmissen&#39;s crypto-attacks](https://github.com/jvdsn/crypto-attacks) and some ideas are inspired from [Yansong Feng](https://github.com/fffmath).

## Requirements

- [**SageMath**](https://www.sagemath.org/) 9.5 with Python 3.10

You can check your SageMath Python version using the following command:

```commandline
$ sage -python --version
Python 3.10.12
```

Note: If your SageMath Python version is older than 3.9.0, some features in given scripts might not work.

## Usage

The standard way to run the attack with the specific parameters requires passing them as command line arguments `sage -python attack.py <modulus_bit_length> <delta> <s>`. For instance, to run the attack with $\ell=512$, $\delta=0.47$, and $s=5$, please run `sage -python attack.py 512 0.47 5`:

```commandline
cPRSA$ sage -python attack.py 512 0.47 5
The parameters:
N = 9247845261559280687272202032856986243220205780672107443336241590071138519185229589720260967480383383456564133610490309706623534189713621588337819964701137
e = 33094679852655913644500990161144486679843079125449690700380916262692468023936108476975980487875653026224553503993694294134942746371316764731541881677400324930953912204129214980874599972506908601054861064166927164487966402091863667740953161390733830061283431939587517936627044851208850457347169039777765518966
Found primes:
113807295748560436061416361496942361650794214297434390995454537370838873516043
81258808591594690780735224428912231883254014636867653291533265583701887722259
The attack costs 3.578 seconds...
```

For instance, to run the attack with $\ell=256$, $\delta=0.5$, and $s=8$, please run `sage -python attack.py 256 0.5 8`:

```commandline
cPRSA$ sage -python attack.py 256 0.5 8
The parameters:
N = 98755123564982306724914430882029374387608510127156928927235354992688559403809
e = 3399965343911515555480161945105754822479061246278335916942324550087062233813130141298408752640694448417369315076174365226930565742698950541880913714076653
Found primes:
339246339697075241646617988151503879727
291101515356552300374416299778317672367
The attack costs 23.718 seconds...
```

## Notes

All the details of the numerical attack experiments are recorded in the `attack.log` file. A more efficient lattice reduction algorithm [flatter](https://github.com/keeganryan/flatter) is used and `USE_FLATTER = True`.

[^cPRSA]: Zheng M., Kunihiro N., Yao Y., "Cryptanalysis of the RSA variant based on cubic Pell equation" | [PDF](https://mengcezheng.github.io/docs/ZKY21.pdf)
