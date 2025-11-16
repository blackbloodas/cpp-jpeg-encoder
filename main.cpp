#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdint.h>
#include <algorithm>
#include <zlib.h>
#include <cmath>
#include <map>
#include <queue>
using namespace std;
uint32_t Width, Height, BytesPerPixel;
const unsigned char pngSignature[] = {0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A};
const unsigned char IHDRSignature[] = {'I', 'H', 'D', 'R'};
const unsigned char IDATSignature[] = {'I', 'D', 'A', 'T'};
const unsigned char IENDSignature[] = {'I', 'E', 'N', 'D'};
const double pi = M_PI;
int Q_Luminance[8][8] = {
    {8, 6, 5, 8, 12, 20, 26, 31},
    {6, 6, 7, 10, 13, 29, 30, 28},
    {7, 7, 8, 12, 20, 29, 35, 28},
    {7, 9, 11, 15, 26, 44, 40, 31},
    {9, 11, 19, 28, 34, 55, 52, 39},
    {12, 18, 28, 32, 41, 52, 57, 46},
    {25, 32, 39, 44, 52, 61, 60, 51},
    {36, 46, 48, 49, 56, 50, 52, 50}};
int Q_Chrominance[8][8] = {
    {9, 9, 12, 24, 50, 50, 50, 50},
    {9, 11, 13, 33, 50, 50, 50, 50},
    {12, 13, 28, 50, 50, 50, 50, 50},
    {24, 33, 50, 50, 50, 50, 50, 50},
    {50, 50, 50, 50, 50, 50, 50, 50},
    {50, 50, 50, 50, 50, 50, 50, 50},
    {50, 50, 50, 50, 50, 50, 50, 50},
    {50, 50, 50, 50, 50, 50, 50, 50}};
// clang-format off
int zigzag[64][2] = {
    {0,0}, {0,1}, {1,0}, {2,0}, {1,1}, {0,2}, {0,3}, {1,2},
    {2,1}, {3,0}, {4,0}, {3,1}, {2,2}, {1,3}, {0,4}, {0,5},
    {1,4}, {2,3}, {3,2}, {4,1}, {5,0}, {6,0}, {5,1}, {4,2},
    {3,3}, {2,4}, {1,5}, {0,6}, {0,7}, {1,6}, {2,5}, {3,4},
    {4,3}, {5,2}, {6,1}, {7,0}, {7,1}, {6,2}, {5,3}, {4,4},
    {3,5}, {2,6}, {1,7}, {2,7}, {3,6}, {4,5}, {5,4}, {6,3},
    {7,2}, {7,3}, {6,4}, {5,5}, {4,6}, {3,7}, {4,7}, {5,6},
    {6,5}, {7,4}, {7,5}, {6,6}, {5,7}, {6,7}, {7,6}, {7,7}
};
// clang-format on
int Bigendianint(int start, const vector<unsigned char> &buffer)
{
    int Value = (static_cast<uint32_t>(buffer[start]) << 24) | // byte << 24 = byte * 2^24 = byte * 16777216
                (static_cast<uint32_t>(buffer[start + 1]) << 16) |
                (static_cast<uint32_t>(buffer[start + 2]) << 8) |
                (static_cast<uint32_t>(buffer[start + 3]));
    return Value;
}
vector<unsigned char> collectIDATData(const vector<unsigned char> &buffer)
{
    int i = 8;
    int chunkLength;
    bool IEND = false;
    vector<unsigned char> Compressed;

    while (!IEND)
    {

        chunkLength = Bigendianint(i, buffer);
        bool IDAT = equal(buffer.begin() + i + 4, buffer.begin() + i + 8, IDATSignature);
        IEND = equal(buffer.begin() + i + 4, buffer.begin() + i + 8, IENDSignature);
        if (IDAT)
        {
            for (int j = 0; j < chunkLength; j++)
            {
                Compressed.push_back(buffer[j + i + 8]);
            }
        }
        i += 12 + chunkLength;
    }

    return Compressed;
}
vector<unsigned char> Deflate(const vector<unsigned char> &Compressed)
{
    uLongf Raw_DataLen = Width * Height * BytesPerPixel + Height;
    vector<unsigned char> Raw_Data(Raw_DataLen);

    int err_uc = uncompress(Raw_Data.data(), &Raw_DataLen, Compressed.data(), Compressed.size());

    if (err_uc == Z_OK)
    {
        Raw_Data.resize(Raw_DataLen);
        return Raw_Data;
    }
    else
    {
        cout << "Decompression failed: " << err_uc << endl;
        return {};
    }
}
int BytesPP(uint8_t ColorType)
{
    switch (ColorType)
    {
    case 0:
        return 1; // Grayscale
    case 2:
        return 3; // RGB
    case 4:
        return 2; // Grayscale + Alpha
    case 6:
        return 4; // RGBA
    default:
        return 0;
    }
}
vector<vector<uint8_t>> UnFilter(const vector<unsigned char> &Filtered)
{
    int Size = Width * BytesPerPixel;
    vector<vector<uint8_t>> Pixels(Height, vector<uint8_t>(Size));

    for (int i = 0; i < Height; i++)
    {
        int Start = i * (Size + 1);
        for (int j = 0; j < Size; j++)
        {
            Pixels[i][j] = Filtered[Start + j + 1];
        }
    }

    for (int i = 0; i < Height; i++)
    {
        if (Filtered[i * (Size + 1)] == 1) // Sub
        {

            for (int j = BytesPerPixel; j < Size; j++)
            {

                Pixels[i][j] += Pixels[i][j - BytesPerPixel];
            }
        }
        else if (Filtered[i * (Size + 1)] == 2 && i > 0) // Up
        {

            for (int j = 0; j < Size; j++)
            {

                Pixels[i][j] += Pixels[i - 1][j];
            }
        }
        else if (Filtered[i * (Size + 1)] == 3) // Average
        {

            for (int j = 0; j < Size; j++)
            {
                if (i == 0 && j < BytesPerPixel)
                {
                    continue;
                }
                else if (i == 0)
                {
                    Pixels[i][j] += Pixels[i][j - BytesPerPixel] / 2;
                }
                else if (j < BytesPerPixel)
                {
                    Pixels[i][j] += Pixels[i - 1][j] / 2;
                }
                else
                {
                    Pixels[i][j] += (Pixels[i][j - BytesPerPixel] + Pixels[i - 1][j]) / 2;
                }
            }
        }
        else if (Filtered[i * (Size + 1)] == 4) // Paeth
        {
            int a, b, c, p;
            int pa, pb, pc;
            int smallestD, smallest;

            for (int j = 0; j < Size; j++)
            {
                a = 0;
                b = 0;
                c = 0;
                if (j >= BytesPerPixel)
                {
                    a = Pixels[i][j - BytesPerPixel];
                }
                if (i > 0)
                {
                    b = Pixels[i - 1][j];
                }
                if (i > 0 && j >= BytesPerPixel)
                {
                    c = Pixels[i - 1][j - BytesPerPixel];
                }
                p = a + b - c;
                pa = abs(p - a);
                pb = abs(p - b);
                pc = abs(p - c);
                smallestD = pb;
                smallest = b;
                if (pa <= pb)
                {
                    smallestD = pa;
                    smallest = a;
                }
                if (pc <= smallestD)
                {
                    smallestD = pc;
                    smallest = c;
                }
                Pixels[i][j] += smallest;
            }
        }
    }
    return Pixels;
}
vector<vector<uint8_t>> YCbCr(const vector<vector<uint8_t>> &Pixels)
{
    vector<vector<uint8_t>> YPixels(Height, vector<uint8_t>(Width * BytesPerPixel));

    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width * BytesPerPixel; j += BytesPerPixel)
        {
            double R = Pixels[i][j];
            double G = Pixels[i][j + 1];
            double B = Pixels[i][j + 2];

            double Y = 0.299 * R + 0.587 * G + 0.114 * B;
            YPixels[i][j] = static_cast<uint8_t>(clamp(Y, 0.0, 255.0));

            double Cb = -0.168736 * R - 0.331264 * G + 0.5 * B + 128.0;
            YPixels[i][j + 1] = static_cast<uint8_t>(clamp(Cb, 0.0, 255.0));

            double Cr = 0.5 * R - 0.418688 * G - 0.081312 * B + 128.0;
            YPixels[i][j + 2] = static_cast<uint8_t>(clamp(Cr, 0.0, 255.0));
        }
    }
    return YPixels;
}
void ChromaDown(vector<vector<uint8_t>> &YPixels)
{

    for (int i = 0; i < Height - 1; i += 2)
    {
        for (int j = 1; j < Width * BytesPerPixel; j += BytesPerPixel * 2)
        {
            double avCb = (YPixels[i][j] + YPixels[i][j + BytesPerPixel] + YPixels[i + 1][j] + YPixels[i + 1][j + BytesPerPixel]) / 4.0;
            YPixels[i][j] = avCb;
            YPixels[i][j + BytesPerPixel] = avCb;
            YPixels[i + 1][j] = avCb;
            YPixels[i + 1][j + BytesPerPixel] = avCb;
            double avCr = (YPixels[i][j + 1] + YPixels[i][j + BytesPerPixel + 1] + YPixels[i + 1][j + 1] + YPixels[i + 1][j + BytesPerPixel + 1]) / 4.0;
            YPixels[i][j + 1] = avCr;
            YPixels[i][j + BytesPerPixel + 1] = avCr;
            YPixels[i + 1][j + 1] = avCr;
            YPixels[i + 1][j + BytesPerPixel + 1] = avCr;
        }
    }
}
vector<vector<vector<uint8_t>>> Blockify(const vector<vector<uint8_t>> &plane)
{
    int blocksHigh = Height / 8;
    int blocksWide = Width / 8;
    int numBlocks = blocksHigh * blocksWide;
    vector<vector<vector<uint8_t>>> Block(numBlocks, vector<vector<uint8_t>>(8, vector<uint8_t>(8)));

    for (int i = 0; i < blocksHigh * 8; i++)
    {
        for (int j = 0; j < blocksWide * 8; j++)
        {
            Block[(i / 8) * blocksWide + (j / 8)][i % 8][j % 8] = plane[i][j];
        }
    }
    return Block;
}
vector<vector<vector<double>>> DCT_Coeff(vector<vector<vector<uint8_t>>> &Block, vector<vector<vector<vector<double>>>> &Basis)
{
    int BlockNum = Block.size();
    vector<vector<vector<double>>> DCT_Coefficient(BlockNum, vector<vector<double>>(8, vector<double>(8, 0.0)));

    for (int m = 0; m < BlockNum; m++)
    {
        for (int u = 0; u < 8; u++)
        {
            for (int v = 0; v < 8; v++)
            {
                for (int i = 0; i < 8; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {

                        DCT_Coefficient[m][u][v] += static_cast<double>(Block[m][i][j]) * Basis[u][v][i][j];
                    }
                }
            }
        }
    }
    return DCT_Coefficient;
}
void Quantization(vector<vector<vector<double>>> &DCT_Block, bool type)
{

    int BlockNum = DCT_Block.size();
    int (*Q)[8];
    if (type)
    {
        Q = Q_Chrominance;
    }
    else
    {
        Q = Q_Luminance;
    }
    for (int m = 0; m < BlockNum; m++)
    {
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < 8; j++)
            {

                DCT_Block[m][i][j] = round(DCT_Block[m][i][j] / Q[i][j]);
            }
        }
    }
}
vector<vector<double>> Ziggify(vector<vector<vector<double>>> &DCT_arr)
{
    int BlockNum = DCT_arr.size();
    vector<vector<double>> zig_arr(BlockNum, vector<double>(64));

    for (int m = 0; m < BlockNum; m++)
    {
        for (int i = 0; i < 64; i++)
        {

            zig_arr[m][i] = DCT_arr[m][zigzag[i][0]][zigzag[i][1]];
        }
    }
    return zig_arr;
}
vector<vector<double>> DPCM(const vector<vector<double>> &zig_arr)
{
    int BlockNum = zig_arr.size();
    vector<vector<double>> DPCM_arr = zig_arr;

    for (int m = 1; m < BlockNum; m++)
    {
        DPCM_arr[m][0] -= zig_arr[m - 1][0];
    }
    return DPCM_arr;
}
vector<vector<int>> RLE(const vector<vector<double>> &DPCM_arr)

{
    int BlockNum = DPCM_arr.size();
    vector<vector<int>> RLE_arr(BlockNum);

    for (int i = 0; i < BlockNum; i++)
    {
        int Zero_N = 0;
        for (int j = 0; j < 64; j++)
        {
            if (DPCM_arr[i][j] == 0)
            {
                Zero_N++;
            }
            else
            {
                RLE_arr[i].push_back(Zero_N);
                RLE_arr[i].push_back(DPCM_arr[i][j]);
                Zero_N = 0;
            }
        }
        RLE_arr[i].push_back(0); // EOB
        RLE_arr[i].push_back(0);
    }

    return RLE_arr;
}
struct Node
{
    pair<int, int> symbol;
    int freq;
    Node *left;
    Node *right;
};

class Comparator
{
public:
    bool operator()(Node *a, Node *b)
    {
        return a->freq > b->freq;
    }
};
static map<pair<int, int>, string> huffmanCodes;
void Traversal(Node *x, string CurrentCode)
{
    string Code = CurrentCode;
    if (x == NULL)
    {
        return;
    }
    else if (x->left == NULL && x->right == NULL)
    {
        huffmanCodes.insert({x->symbol, Code});
        ;
        return;
    }
    else
    {
        Traversal(x->left, CurrentCode + "0");
        Traversal(x->right, CurrentCode + "1");
    }
}
void Huffman(vector<vector<int>> &RLE)
{
    int BlockNum = RLE.size();
    map<pair<int, int>, int> Frequency;

    for (int i = 0; i < BlockNum; i++)
    {
        for (int j = 0; j < RLE[i].size(); j += 2)
        {
            int zero_run = RLE[i][j];
            int value = RLE[i][j + 1];
            Frequency[{zero_run, value}]++;
        }
    }

    priority_queue<Node *, vector<Node *>, Comparator> pq;
    for (const auto &pair : Frequency)
    {
        Node *node = new Node();
        node->symbol = pair.first;
        node->freq = pair.second;
        pq.push(node);
    }
    while (pq.size() > 1)
    {
        Node *NodeParent = new Node();
        Node *Leaf1 = pq.top();
        pq.pop();
        Node *Leaf2 = pq.top();
        pq.pop();
        NodeParent->freq = Leaf1->freq + Leaf2->freq;
        NodeParent->right = Leaf1;
        NodeParent->left = Leaf2;
        pq.push(NodeParent);
    }
    Node *root = pq.top();
    Traversal(root, "");
}
vector<unsigned char> EncodeRLE(const vector<vector<int>> &RLE)
{
    int BlockNum = RLE.size();
    vector<unsigned char> HRLE;

    string bitString;

    for (int i = 0; i < BlockNum; i++)
    {
        for (int j = 0; j < RLE[i].size(); j += 2)
        {
            int Zero_Run = RLE[i][j];
            int Value = RLE[i][j + 1];
            bitString += huffmanCodes[{Zero_Run, Value}];
        }
    }
    for (int i = 0; i < bitString.length() / 8; i++)
    {
        string binaryString = bitString.substr(i * 8, 8);
        unsigned char integerValue = stoi(binaryString, nullptr, 2);
        HRLE.push_back(integerValue);
    }

    int remaining = bitString.length() % 8;
    if (remaining > 0)
    {
        string lastBits = bitString.substr(bitString.length() - remaining);
        lastBits.append(8 - remaining, '0');
        unsigned char byteValue = stoi(lastBits, nullptr, 2);
        HRLE.push_back(byteValue);
    }
    return HRLE;
}

void DisplayASCII(const vector<vector<uint8_t>> &Pixels)
{
    const char chars[] = " .:-=+*#%@";

    for (int i = 0; i < Height; i += 16)
    {
        for (int j = 0; j < Width * BytesPerPixel; j += BytesPerPixel * 6)
        {
            int gray = (Pixels[i][j] + Pixels[i][j + 1] + Pixels[i][j + 2]) / 3;
            cout << chars[gray * 10 / 256];
        }
        cout << endl;
    }
}

int main()
{

    string signiture;
    const string Image = "Ikebukuro-Cat-Park-3.png";
    ifstream MyFile(Image, ios::in | ios::binary);
    if (!MyFile)
    {
        cout << "Failed to open file!" << endl;
        return 1;
    }
    else
    {
        cout << "File loaded" << endl;
    }

    vector<unsigned char> buffer(istreambuf_iterator<char>(MyFile), {});
    MyFile.close();

    bool validSignature = true;

    for (int i = 0; i < 8; ++i)
    {
        if (buffer[i] != pngSignature[i])
        {
            validSignature = false;
            break;
        }
    }

    if (!validSignature)
    {
        cout << "Invalid File Type" << endl;
    }
    else
    {
        cout << "PNG: " << endl;

        for (size_t i = 0; i < 16; ++i)
        {

            cout << hex << setw(2) << setfill('0') << static_cast<uint32_t>(buffer[i]) << " ";
            if ((i + 1) % 4 == 0)
            {
                cout << endl;
            }
        }

        uint8_t BitDepth = buffer[24];
        uint8_t ColorType = buffer[25];
        uint8_t InterlaceMethod = buffer[28];
        Width = Bigendianint(16, buffer);
        Height = Bigendianint(20, buffer);
        BytesPerPixel = BytesPP(ColorType);

        cout << dec;
        cout << "Width: " << Width << endl;
        cout << "Height: " << Height << endl;
        cout << "Bit Depth: " << static_cast<int>(BitDepth) << endl;
        cout << "Color Type: " << static_cast<int>(ColorType) << endl;
        cout << "Interlace Method: " << static_cast<int>(InterlaceMethod) << endl;
        vector<unsigned char> Compressed = collectIDATData(buffer);
        cout << "Collected " << Compressed.size() / 1024 << "Kbytes of compressed data" << endl;
        vector<unsigned char> Raw_Data = Deflate(Compressed);
        if (!Raw_Data.empty())
        {
            cout << "Decompressed " << Raw_Data.size() / (1024 * 1024) << "MB of raw image data" << endl;
        }

        vector<vector<uint8_t>> Pixels = UnFilter(Raw_Data);
        DisplayASCII(Pixels);
        vector<vector<uint8_t>> YPixels = YCbCr(Pixels);
        ChromaDown(YPixels);
        vector<vector<uint8_t>> Y_plane(Height, vector<uint8_t>(Width));
        vector<vector<uint8_t>> Cb_plane(Height, vector<uint8_t>(Width));
        vector<vector<uint8_t>> Cr_plane(Height, vector<uint8_t>(Width));
        for (int i = 0; i < Height; i++)
        {
            for (int j = 0; j < Width; j++)
            {
                Y_plane[i][j] = YPixels[i][j * BytesPerPixel];
                Cb_plane[i][j] = YPixels[i][(j * BytesPerPixel) + 1];
                Cr_plane[i][j] = YPixels[i][(j * BytesPerPixel) + 2];
            }
        }
        vector<vector<vector<uint8_t>>> Y_blocks = Blockify(Y_plane);
        vector<vector<vector<uint8_t>>> Cb_blocks = Blockify(Cb_plane);
        vector<vector<vector<uint8_t>>> Cr_blocks = Blockify(Cr_plane);
        vector<vector<vector<vector<double>>>> Basis(8, vector<vector<vector<double>>>(8, vector<vector<double>>(8, vector<double>(8))));

        for (int u = 0; u < 8; u++)
        {
            for (int v = 0; v < 8; v++)
            {
                for (int i = 0; i < 8; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        double a_u = 1 / sqrt(8);
                        double a_v = 1 / sqrt(8);

                        if (u > 0)
                        {
                            a_u = sqrt(2.0) / sqrt(8.0);
                        }
                        if (v > 0)
                        {
                            a_v = sqrt(2.0) / sqrt(8.0);
                        }
                        Basis[u][v][i][j] = a_u * a_v * cos(((2 * i + 1) * u * pi) / 16.0) * cos(((2 * j + 1) * v * pi) / 16.0);
                    }
                }
            }
        }
        vector<vector<vector<double>>> Y_DCT = DCT_Coeff(Y_blocks, Basis);
        vector<vector<vector<double>>> Cb_DCT = DCT_Coeff(Cb_blocks, Basis);
        vector<vector<vector<double>>> Cr_DCT = DCT_Coeff(Cr_blocks, Basis);
        Quantization(Y_DCT, 0);
        Quantization(Cb_DCT, 1);
        Quantization(Cr_DCT, 1);
        vector<vector<double>> Y_Zigzag = Ziggify(Y_DCT);
        vector<vector<double>> Cb_Zigzag = Ziggify(Cb_DCT);
        vector<vector<double>> Cr_Zigzag = Ziggify(Cr_DCT);
        vector<vector<double>> Y_DCPM = DPCM(Y_Zigzag);
        vector<vector<double>> Cb_DCPM = DPCM(Cb_Zigzag);
        vector<vector<double>> Cr_DCPM = DPCM(Cr_Zigzag);
        vector<vector<int>> Y_RLE = RLE(Y_DCPM);
        vector<vector<int>> Cb_RLE = RLE(Cb_DCPM);
        vector<vector<int>> Cr_RLE = RLE(Cr_DCPM);

        int original_size = Height * Width * BytesPerPixel;
        int rle_size = 0;
        for (const auto &block : Y_RLE)
            rle_size += block.size();
        for (const auto &block : Cb_RLE)
            rle_size += block.size();
        for (const auto &block : Cr_RLE)
            rle_size += block.size();

        cout << "\n=== Compression Statistics ===" << endl;
        cout << "Original size: " << original_size / 1024 << " KB" << endl;
        cout << "RLE encoded size: " << rle_size * sizeof(int) / 1024 << " KB" << endl;
        cout << "Compression ratio (): " << (float)original_size / (rle_size * sizeof(int)) << ":1" << endl;

        Huffman(Y_RLE);
        vector<unsigned char> Y_Huffman = EncodeRLE(Y_RLE);

        huffmanCodes.clear();
        Huffman(Cb_RLE);
        vector<unsigned char> Cb_Huffman = EncodeRLE(Cb_RLE);

        huffmanCodes.clear();
        Huffman(Cr_RLE);
        vector<unsigned char> Cr_Huffman = EncodeRLE(Cr_RLE);

        int huffman_size = Y_Huffman.size() + Cb_Huffman.size() + Cr_Huffman.size();

        cout << "\n=== Final Huffman Compression ===" << endl;
        cout << "Y channel Huffman: " << Y_Huffman.size() << " bytes" << endl;
        cout << "Cb channel Huffman: " << Cb_Huffman.size() << " bytes" << endl;
        cout << "Cr channel Huffman: " << Cr_Huffman.size() << " bytes" << endl;
        cout << "Total Huffman encoded size: " << huffman_size / 1024 << " KB" << endl;
        cout << "Final compression ratio: " << (float)original_size / huffman_size << ":1" << endl;
        cout << "Space saved: " << (1.0 - (float)huffman_size / original_size) * 100 << "%" << endl;
    }

    return 0;
}
