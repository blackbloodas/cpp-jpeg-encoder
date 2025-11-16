Implementing JPEG Compression from Scratch in C++

This project implements a complete JPEG-style image compression pipeline in C++ from scratch, without relying on image processing libraries. Starting from a PNG file, the system decodes the image data, applies a series of signal processing and compression algorithms, and achieves a compression ratio of 31.7:1 with 96.8% space savings. Every component from PNG chunk parsing to Huffman encoding was implemented at the algorithmic level to deeply understand the mathematics and data structures behind modern image compression.
