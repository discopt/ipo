----Installation IPO Wrapper----

Requirements:
1) Python (ab 2.7) and Cython (0.24.1) installed
2) ScipOptSuite installed in /opt
3) IPO installed as shared library
4) Sage is installed

Installation:
1)
>>python setup.py build
in this folder

2)
/build/lib.linux-x86_64-2.7/IPO.so
copy to Sage root folder (SageMath folder)

Start:
1) ./sage
2) import IPO
