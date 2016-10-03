#!/bin/bash
cd modeling/DAECalc/daskr/src
# f2py -c daskr.pyf ddaskr.f daux.f dlinpk.f# mv _daskr.so ..


cd ../../../..
f2py -c __PC12_48__.pyf __PC12_48__.c
f2py -c __dPC12_48__.pyf __dPC12_48__.c
f2py -c __d2PC12_48__.pyf __d2PC12_48__.c
f2py -c __PC12_MA__.pyf __PC12_MA__.c 
f2py -c __dPC12_MA__.pyf __dPC12_MA__.c 
f2py -c __d2PC12_MA__.pyf __d2PC12_MA__.c 
