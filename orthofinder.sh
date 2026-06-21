#!/bin/bash

orthofinder -f data/genomes -o results

#
# Orthofinder feels very straightforward.
#
# I installed with miniconda (conda 25.9.1, python 3.13),
# by making an environment while requesting for these packages:
# # conda create -n ortho -c bioconda orthofinder famsa mcl fastme numpy scipy scikit-learn biopython
# # conda activate ortho
# 
# Apparently, that downloaded these packages in the following versions:
# # Name                     Version          Build             Channel
# _libgcc_mutex              0.1              main
# _openmp_mutex              5.1              1_gnu
# aria2                      1.23.0           0                 bioconda
# biopython                  1.70             np112py35_1       bioconda
# blas                       1.0              mkl
# blast                      2.5.0            hc0b0e79_3        bioconda
# boost                      1.67.0           py35_4
# bzip2                      1.0.8            h5eee18b_6
# c-ares                     1.11.0           h470a237_1        bioconda
# ca-certificates            2025.12.2        h06a4308_0
# certifi                    2020.6.20        pyhd3eb1b0_3
# diamond                    2.1.11           h5ca1c30_2        bioconda
# famsa                      2.4.1            h9ee0642_0        bioconda
# fastme                     2.1.6.3          h7b50bb2_1        bioconda
# fasttree                   2.2.0            h7b50bb2_1        bioconda
# freetype                   2.14.1           hf5b9546_0
# gawk                       5.3.1            h5eee18b_0
# gmp                        6.3.0            h6a678d5_0
# icu                        58.2             he6710b0_3
# intel-openmp               2025.0.0         h06a4308_1171
# iqtree                     3.0.1            h503566f_0        bioconda
# jpeg                       9f               h5ce9db8_0
# lerc                       4.0.0            h6a678d5_0
# libboost                   1.67.0           h46d08c1_4
# libdeflate                 1.22             h5eee18b_0
# libffi                     3.3              he6710b0_2
# libgcc                     15.2.0           h69a1729_7
# libgcc-ng                  15.2.0           h166f726_7
# libgfortran-ng             7.5.0            ha8ba4b0_17
# libgfortran4               7.5.0            ha8ba4b0_17
# libgomp                    15.2.0           h4751f2c_7
# libnsl                     2.0.0            h5eee18b_0
# libpng                     1.6.56           h22898a0_0
# libssh2                    1.10.0           h37d81fd_2
# libstdcxx                  15.2.0           h39759b7_7
# libstdcxx-ng               15.2.0           hc03a8fd_7
# libtiff                    4.5.1            hffd6297_1
# libwebp-base               1.6.0            hb7bb969_0
# libxcb                     1.17.0           h9b100fa_0
# libxml2                    2.10.4           hcbfbd50_0
# libzlib                    1.3.1            hb25bd0a_0
# lz4-c                      1.9.4            h6a678d5_1
# mafft                      7.525            h031d066_1        bioconda
# mcl                        22.282           pl5321h7b50bb2_4  bioconda
# mkl                        2018.0.3         1
# mmseqs2                    18.8cc5c         hd6d6fdc_0        bioconda
# mmtf-python                1.0.2            py35_0            bioconda
# mpi                        1.0              openmpi
# msgpack-python             0.5.6            py35h6bb024c_1
# muscle                     5.3              h9948957_3        bioconda
# ncurses                    6.5              h7934f7d_0
# numpy                      1.12.1           py35hca0bb5e_1
# olefile                    0.46             pyhd3eb1b0_0
# openmpi                    4.0.2            hb1b8bf9_1
# openssl                    1.1.1w           h7f8727e_0
# orthofinder                2.5.4            hdfd78af_0        bioconda
# perl                       5.32.1           0_h5eee18b_perl5
# pillow                     5.2.0            py35heded4f4_0
# pip                        10.0.1           py35_0
# pthread-stubs              0.3              h0ce48e5_1
# py-boost                   1.67.0           py35h04863e7_4
# python                     3.5.6            h12debd9_1
# raxml                      8.2.13           h7b50bb2_3        bioconda
# raxml-ng                   0.9.0            h192cbe9_1        bioconda
# readline                   8.3              hc2a1206_0
# reportlab                  3.5.6            py35he686d34_0
# scikit-learn               0.20.0           py35h4989274_1
# scipy                      1.1.0            py35hd20e5f9_0
# setuptools                 40.2.0           py35_0
# sqlite                     3.51.2           h3e8d24a_0
# tk                         8.6.15           h54e0aa7_0
# wheel                      0.37.1           pyhd3eb1b0_0
# xorg-libx11                1.8.12           h9b100fa_1
# xorg-libxau                1.0.12           h9b100fa_0
# xorg-libxdmcp              1.1.5            h9b100fa_0
# xorg-xorgproto             2024.1           h5eee18b_1
# xz                         5.8.2            h448239c_0
# zlib                       1.3.1            hb25bd0a_0
# zstd                       1.5.7            h11fc155_0


