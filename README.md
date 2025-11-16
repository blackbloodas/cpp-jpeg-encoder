Implementing JPEG Compression from Scratch in C++

Abdelrahman Elkaseer

1. Brief
This project implements a complete JPEG-style image compression pipeline in C++ from scratch, without relying on image processing libraries. Starting from a PNG file, the system decodes the image data, applies a series of signal processing and compression algorithms, and achieves a compression ratio of 31.7:1 with 96.8% space savings. Every component from PNG chunk parsing to Huffman encoding was implemented at the algorithmic level to deeply understand the mathematics and data structures behind modern image compression.

2. Pipeline


2.1. PNG Decoding
1. PNG Signature & Chunk Identifiers
PNG files use specific byte sequences to identify file structure. The 8-byte signature validates the file format, while 4-byte codes identify different chunk types (IHDR for headers, IDAT for image data, IEND for file end).
2. Big-Endian Integer Conversion
PNG stores multi-byte integers in big-endian format (most significant byte first). This function reconstructs 32-bit integers by shifting each byte to its proper position and combining them with bitwise OR operations.
3. IDAT Chunk Collection
This function iterates through the PNG chunk structure, identifying and extracting all IDAT chunks that contain the compressed pixel data. Each chunk has a length field, type identifier, data payload, and CRC checksum. The loop continues until reaching the IEND marker
4. DEFLATE Decompression
PNG uses DEFLATE compression (same algorithm as ZIP). The zlib library decompresses the collected IDAT data. The buffer size accounts for all pixels plus one filter-type byte per scanline.
5. Color Type Handler
PNG supports different color formats with varying bytes per pixel. This function maps color type codes to their corresponding byte counts essential for correctly interpreting the pixel data buffer.


2.2. Unfilter
Adaptive Filter Reversal
PNG applies adaptive filtering to each scanline before compression to improve DEFLATE efficiency. Each scanline begins with a filter type byte (0-4) that indicates which prediction algorithm was used. The unfiltering process reverses these predictions by adding back the predicted values to reconstruct the original pixels.
Filter Types Implemented:
Filter 0 (None): No filtering applied pixels stored as-is.
Filter 1 (Sub): Each pixel is the difference from the pixel to its left. Reversal adds the left neighbor.
Filter 2 (Up): Each pixel is the difference from the pixel above it. Reversal adds the upper neighbor.
Filter 3 (Average): Each pixel is the difference from the average of left and upper neighbors. Reversal adds the floor of their average.
Filter 4 (Paeth): Uses the Paeth predictor a non-linear function that selects the neighboring pixel (left, up, or upper-left) closest to a predicted value. The predictor computes p = a + b - c and chooses whichever of a, b, or c is nearest to p by absolute distance.



2.3. Color Space Transformation
RGB to YCbCr Conversion
JPEG compression exploits human visual perception by separating luminance (brightness) from chrominance (color). The YCbCr color space represents images as one luminance channel (Y) and two chrominance channels (Cb and Cr), allowing independent processing and more aggressive compression of color information.

Why This Matters:
Human vision is more sensitive to brightness variations than color changes. By separating Y from Cb/Cr, we can downsample the chrominance channels without perceptually degrading the image. This is the first step toward lossy compression focusing precision where human eyes are most sensitive.


2.4. Chroma Subsampling
4:2:0 Chroma Downsampling
Chroma subsampling exploits human visual system limitations our eyes are far more sensitive to brightness (luminance) than color (chrominance). This stage reduces chrominance resolution while preserving full luminance resolution, achieving significant data reduction with minimal perceptual quality loss.

4:2:0 Subsampling Scheme:
The implementation averages chrominance values over 2×2 pixel blocks:
	Luminance (Y): Full resolution maintained—every pixel keeps its brightness value
	Cb and Cr channels: Reduced to 1/4 resolution by averaging four neighboring pixels
	Each 2×2 block shares a single Cb value and single Cr value
________________________________________

Why 4:2:0?
The notation describes chrominance sampling relative to luminance in a 4×2 pixel region:
	4: Four luminance samples per row (full horizontal resolution)
	2: Two chrominance samples per row (half horizontal resolution)
	0: No additional chrominance samples in the second row (half vertical resolution)
	Data Reduction:
	This step reduces color data by 75% with negligible perceptual impact. The human eye's reduced sensitivity to color detail makes this lossy compression nearly invisible in most images.

	
2.5. Block Partitioning
8×8 Block Decomposition
The Discrete Cosine Transform (DCT) operates on small blocks rather than entire images. This function partitions each color plane into non-overlapping 8×8 pixel blocks the fundamental unit for JPEG compression. This block size balances compression efficiency with computational complexity.

Partitioning Process:
	Calculate block dimensions: blocksWide = Width / 8, blocksHigh = Height / 8
	Total blocks: blocksWide × blocksHigh per channel
	Raster-to-block mapping: Pixels are reorganized from a 2D image array into a 3D structure where each element is an 8×8 block
________________________________________

Why 8×8 Blocks?
	Computational efficiency: Small enough for fast DCT computation
	Locality: Captures local image features and correlations
	Standard compliance: JPEG baseline specification requires 8×8 blocks
	Practical balance: Larger blocks would capture more context but increase complexity; smaller blocks would reduce compression efficiency


2.6. Discrete Cosine Transform (DCT)
Frequency Domain Transformation
The DCT is the mathematical heart of JPEG compression. It transforms spatial pixel values into frequency coefficients, concentrating image energy into a few low-frequency components while spreading noise and detail across high-frequency components that can be aggressively quantized.
 

DCT Basis Functions:
The transformation decomposes each 8×8 block into a weighted sum of 64 basis patterns. Each basis function represents a specific spatial frequency:
Basis equation:
Basis[u][v][i][j] = α(u) · α(v) · cos[π·u/8·(i + 0.5)] · cos[π·v/8·(j + 0.5)]
Where α(k) = 1/√8 for k=0, and α(k) = √(2/8) otherwise.
	(u,v) = (0,0): DC coefficient—represents the average block intensity
	Low (u,v): Smooth gradients and large-scale features
	High (u,v): Sharp edges, texture, and fine detail

________________________________________
Transformation Process:
For each 8×8 block, every DCT coefficient is computed by multiplying the block pixels by the corresponding basis function and summing the results. This produces 64 frequency coefficients where most image energy concentrates in the upper-left (low frequencies).
________________________________________

Why DCT Works:
Natural images exhibit strong spatial correlation neighboring pixels tend to have similar values. The DCT exploits this by packing energy into low-frequency coefficients, allowing high-frequency coefficients (often near zero) to be heavily compressed with minimal visual impact.


2.7. Quantization
Lossy Compression Step
Quantization is where JPEG becomes lossy. DCT coefficients are divided by quantization table values and rounded to integers, discarding precision based on perceptual importance. This is the primary quality control mechanism larger quantization values mean more aggressive compression and lower quality.
Quantization Tables:
Two standard tables control compression for different channels

Luminance Table (Q_Luminance): 
Used for the Y channel
	Preserves low frequencies (smaller divisors in upper-left)
	Aggressively quantizes high frequencies (larger divisors in lower-right)
	Reflects human sensitivity to brightness detail

Chrominance Table (Q_Chrominance): 
Used for Cb and Cr channels
	More aggressive across all frequencies
	Exploits reduced human sensitivity to color detail
	Higher compression ratios with minimal perceptual loss

Mathematical Operation:
For each DCT coefficient: Quantized[i][j] = round(DCT[i][j] / Q[i][j])
________________________________________
Impact:
High-frequency coefficients often quantize to zero, creating long runs of zeros that RLE will exploit. This step achieves the bulk of JPEG's compression while determining final image quality.

2.8. Zigzag Ordering
Frequency-Sequential Arrangement
Zigzag ordering reorders the 8×8 coefficient matrix into a 1D sequence, progressing from low to high frequencies. This arrangement clusters non-zero coefficients at the beginning and zeros at the end, maximizing run-length encoding efficiency in the next stage.
 

Why Zigzag?
After quantization, high-frequency coefficients (lower-right of the block) are predominantly zero. The zigzag pattern ensures these zeros appear consecutively in the 1D array, creating long runs of zeros that run-length encoding can compress efficiently.


2.9. Differential Pulse Code Modulation (DPCM)
DC Coefficient Encoding
DPCM exploits correlation between adjacent blocks by encoding DC coefficients (position 0 in each zigzag array) as differences rather than absolute values. Since neighboring blocks often have similar average intensities, the differences are typically small numbers that compress more efficiently.

Differential Encoding:
Instead of storing DC values directly, each block stores the difference from the previous block's DC coefficient:
	Block 0: DC value stored as-is (no previous block to reference)
	Block 1: DC_diff = DC[1] - DC[0]
	Block 2: DC_diff = DC[2] - DC[1]
	Block n: DC_diff = DC[n] - DC[n-1]
________________________________________

Why Only DC Coefficients?
The DC coefficient represents the average block intensity and exhibits strong spatial correlation adjacent image regions tend to have similar brightness. AC coefficients (positions 1-63) show less inter-block correlation and are left unchanged.
________________________________________

Compression Benefit:
Differences cluster around zero with smaller magnitude than absolute values, requiring fewer bits to encode. This differential representation is particularly effective in smooth image regions where DC values change gradually.



2.10. Run-Length Encoding (RLE)
Zero Run Compression
Run-length encoding exploits the long sequences of zeros created by quantization and zigzag ordering. Instead of storing each zero individually, RLE represents consecutive zeros as (zero_count, value) pairs, dramatically reducing data volume.




Encoding Format:
Each non-zero coefficient is encoded as a pair:
	First element: Number of zeros preceding this coefficient
	Second element: The coefficient value itself
________________________________________
End-of-Block (EOB) Marker:
The special pair (0, 0) signals that all remaining coefficients in the block are zero, avoiding the need to encode long trailing zero runs. This is crucial since most blocks end with 40-50+ consecutive zeros after quantization.
________________________________________

Compression Achievement:
A 64-element array often compresses to 6-12 pairs, reducing data by 80-90%. High-frequency regions with more detail produce longer RLE sequences, while smooth regions compress to just a few pairs.


2.11. Huffman Encoding
Optimal Entropy Encoding
Huffman coding assigns variable-length binary codes to RLE symbol pairs based on frequency common symbols get short codes, rare symbols get longer codes. This is the final lossless compression stage, creating optimal prefix-free codes that minimize total bit usage.


Algorithm Steps:
1. Frequency Analysis: Count occurrences of each (zero_run, value) pair across all blocks
2. Priority Queue Construction: Create leaf nodes for each unique symbol, prioritized by frequency (min-heap)
3. Tree Building:
	Extract two lowest-frequency nodes
	Create parent node with combined frequency
	Repeat until one root node remains
4. Code Generation: Traverse tree recursively left branches add '0', right branches add '1'. Each leaf's path becomes its binary code.
________________________________________
Prefix-Free Property:
No code is a prefix of another, enabling unambiguous decoding. Reading the bitstream left-to-right, each complete code uniquely identifies one symbol.

2.12. Bit Packing
Final Compression Stage
This function applies the Huffman codes to encode the RLE data and packs the variable-length bit sequences into bytes. The output is a compact binary stream representing the fully compressed image data. 

Encoding Process:
1. Code Lookup & Concatenation: For each (zero_run, value) pair, look up its Huffman code and append to a continuous bit string.
2. Bit Packing: Group the bit string into 8-bit chunks and convert each to a byte:
	"11001101" → byte value 205
	"10000000" → byte value 128 (last bits padded with zeros)
3. Partial Byte Handling: If the total bit length isn't divisible by 8, pad the final byte with zeros to complete it.
________________________________________

Why Bit-Level Encoding?
Huffman codes don't align to byte boundaries—a code might be 3 bits while another is 7 bits. Efficient storage requires packing these variable-length codes tightly into bytes, wasting no space.


3. Results
3.1 Console Output


Performance Metrics
Achieved Compression:
	Compression Ratio: 31.7:1
	Space Savings: 96.8%
	Data Reduction: From megabytes to kilobytes
Stage-by-Stage Contribution:
	Chroma Subsampling: 75% reduction in color data
	DCT + Quantization: Concentrates energy, enables aggressive quantization
	Zigzag + RLE: 80-90% reduction through zero-run compression
	Huffman Coding: Final 20-30% improvement through optimal symbol encoding
