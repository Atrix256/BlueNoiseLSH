// Configurable settings
#define DIMENSION()    2      // dimensionality of the data points
#define NUMPOINTS()    1000   // the number of randomly generated (white noise) data points
#define HASHCOUNT()    8      // how many hashes are done for each point
#define MINHASHCOUNT() 3      // the minimum number of hash collisions that must occur for a point to be considered at all close
#define POINTDOMAIN()  5      // the coordinates go from - this value to + this value

#define IMAGESIZE()       500 // the size of the image - width and height both.

#define DOANGLEIMAGETEST() 1

// TODO: tune HASHCOUNT() and the other parameters!

#include <random>
#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <stdint.h>
#include <complex>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_MSC_SECURE_CRT
#include "stb_image_write.h"

// -------------------------------------------------------------------------------

typedef std::array<float, DIMENSION()> TPoint;
typedef uint32_t uint32;
typedef uint8_t uint8;

enum class PointID : uint32 {};

static const float c_pi = 3.14159265359f;
static const float c_goldenRatioConjugate = 0.61803398875f;  // 1 / golden ratio. golden ratio is the positive solution to x^2 - x - 1 = 0

static const float c_goldenRatio2 = 1.32471795724474602596f;  // positive solution to x^3-x-1 = 0
static const float c_goldenRatio2Conjugate = 1.0f / c_goldenRatio2;

// -------------------------------------------------------------------------------

struct Image
{
    Image(int width, int height)
    {
        m_width = width;
        m_height = height;
        m_pixels.resize(m_width*m_height * 4); // 4 channels per pixel
        std::fill(m_pixels.begin(), m_pixels.end(), 255);
    }

    int m_width;
    int m_height;
    std::vector<uint8> m_pixels;
};

// -------------------------------------------------------------------------------

struct ImageComplex
{
    ImageComplex(int width, int height)
    {
        m_width = width;
        m_height = height;
        m_pixels.resize(m_width*m_height);
        std::fill(m_pixels.begin(), m_pixels.end(), std::complex<float>(0.0f, 0.0f));
    }

    int m_width;
    int m_height;
    std::vector<std::complex<float>> m_pixels;
};

// -------------------------------------------------------------------------------

std::complex<float> DFTPixel (const Image &image, int K, int L)
{
    std::complex<float> ret(0.0f, 0.0f);

    const uint8* pixel = image.m_pixels.data();
    for (int y = 0; y < image.m_height; ++y)
    {
        for (int x = 0; x < image.m_width; ++x)
        {
            float grey = float(pixel[0]) / 255.0f;
            float v = float(K * x) / float(image.m_width);
            v += float(L * y) / float(image.m_height);
            ret += std::complex<float>(grey, 0.0f) * std::polar<float>(1.0f, -2.0f * c_pi * v);
            pixel += 4;
        }
    }

    return ret;
}

void DFTImage (const Image &srcImage, ImageComplex &destImage)
{
    // calculate 2d dft (brute force, not using fast fourier transform)
    destImage = ImageComplex(srcImage.m_width, srcImage.m_height);
    std::complex<float>* pixel = destImage.m_pixels.data();
    for (int y = 0; y < srcImage.m_height; ++y)
    {
        for (int x = 0; x < srcImage.m_width; ++x)
        {
            *pixel = DFTPixel(srcImage, x, y);
            pixel++;
        }
    }
}

void GetMagnitudeData (const ImageComplex& srcImage, Image& destImage)
{
    // size the output image
    destImage = Image(srcImage.m_width, srcImage.m_height);

    // get floating point magnitude data
    std::vector<float> magArray;
    magArray.resize(srcImage.m_width*srcImage.m_height);
    float maxmag = 0.0f;
    for (int x = 0; x < srcImage.m_width; ++x)
    {
        for (int y = 0; y < srcImage.m_height; ++y)
        {
            // Offset the information by half width & height in the positive direction.
            // This makes frequency 0 (DC) be at the image origin, like most diagrams show it.
            int k = (x + srcImage.m_width / 2) % srcImage.m_width;
            int l = (y + srcImage.m_height / 2) % srcImage.m_height;
            const std::complex<float> &src = srcImage.m_pixels[l*srcImage.m_width + k];

            float mag = std::abs(src);
            if (mag > maxmag)
                maxmag = mag;

            magArray[y*srcImage.m_width + x] = mag;
        }
    }
    if (maxmag == 0.0f)
        maxmag = 1.0f;

    const float c = 255.0f / log(1.0f+maxmag);

    // normalize the magnitude data and send it back in [0, 255]
    for (int x = 0; x < srcImage.m_width; ++x)
    {
        for (int y = 0; y < srcImage.m_height; ++y)
        {
            float src = c * log(1.0f + magArray[y*srcImage.m_width + x]);

            uint8 magu8 = uint8(src);

            uint8* dest = &destImage.m_pixels[(y*destImage.m_width + x) * 4];
            dest[0] = magu8;
            dest[1] = magu8;
            dest[2] = magu8;
            dest[3] = 255;
        }
    }
}

// -------------------------------------------------------------------------------

TPoint Multiply(const TPoint& point, const std::array<float, DIMENSION()*DIMENSION()>& matrix)
{
    TPoint ret;
    for (int i = 0; i < DIMENSION(); ++i)
    {
        int rowOffset = i * DIMENSION();

        ret[i] = 0.0f;
        for (int j = 0; j < DIMENSION(); ++j)
            ret[i] += point[j] * matrix[rowOffset + j];
    }
    return ret;
}

// -------------------------------------------------------------------------------

std::array<float, DIMENSION()*DIMENSION()> Transpose(const std::array<float, DIMENSION()*DIMENSION()>& matrix)
{
    std::array<float, DIMENSION()*DIMENSION()> ret;

    for (int i = 0; i < DIMENSION(); ++i)
        for (int j = 0; j < DIMENSION(); ++j)
            ret[i * DIMENSION() + j] = matrix[j * DIMENSION() + i];

    return ret;
}

// -------------------------------------------------------------------------------

struct PointHashData
{
    std::array<float, DIMENSION()*DIMENSION()> rotation;
    float offsetX;
    float angle; // what made the rotation matrix

    std::unordered_map<int, std::vector<PointID>> m_hashBuckets;
};

typedef void(*GeneratePointHashDatas) (std::array<PointHashData, HASHCOUNT()>& hashDatas);

class LHS
{
public:
    LHS (const std::vector<TPoint>& points, GeneratePointHashDatas generatePointHashDatas)
    {
        generatePointHashDatas(m_pointHashDatas);

        m_points = points;
        for (size_t i = 0; i < m_points.size(); ++i)
            AddPoint(m_points[i], PointID(i));
    }

    void Query (const TPoint& point, std::unordered_map<PointID, int>& results) const
    {
        for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
        {
            const PointHashData& pointHashData = m_pointHashDatas[hashIndex];

            int hashValue = HashPoint(point, hashIndex);

            auto it = pointHashData.m_hashBuckets.find(hashValue);

            if (it == pointHashData.m_hashBuckets.end())
                continue;

            for (PointID pointID : it->second)
            {
                auto it = results.find(pointID);
                if (it != results.end())
                {
                    it->second++;
                }
                else
                {
                    results[pointID] = 1;
                }
            }
        }
    }

    const TPoint& GetPoint(PointID point) const
    {
        return m_points[(uint32)point];
    }

    const PointHashData& GetPointHashData (int hashIndex) const
    {
        return m_pointHashDatas[hashIndex];
    }

private:

    int HashPoint (const TPoint& point, int hashIndex) const
    {
        // calculate the hash of the point by rotating it, adding to it's x component and flooring that x component.
        const PointHashData& pointHashData = m_pointHashDatas[hashIndex];
        TPoint rotatedPoint = Multiply(point, pointHashData.rotation);
        return (int)std::floorf(rotatedPoint[0] + pointHashData.offsetX);
    }

    void AddPoint (const TPoint& point, PointID pointID)
    {
        for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
        {
            // store the point in the hash bucket
            int hashValue = HashPoint(point, hashIndex);
            m_pointHashDatas[hashIndex].m_hashBuckets[hashValue].push_back(pointID);
        }
    }

private:

    std::vector<TPoint> m_points; // a PointID (uint32) is the point's index in this list

    std::array<PointHashData, HASHCOUNT()> m_pointHashDatas;
};

// -------------------------------------------------------------------------------

void printfPoint(const char* name, const TPoint& point)
{
    if (name)
        printf("%s: ", name);

    bool first = true;
    for (float f : point)
    {
        if (first)
            printf("(");
        else
            printf(", ");
        printf("%0.2f", f);
        first = false;
    }
    printf(")\n");
}

// -------------------------------------------------------------------------------

std::mt19937& RNG()
{
    static std::random_device rd;
    static std::mt19937 rng(rd());
    return rng;
}

// -------------------------------------------------------------------------------

float SimilarityScore(const PointHashData& A, const PointHashData& B)
{
    // TODO: need to factor in offset somehow!

    // the similarity score of two rotation matrices is the Frobenius inner product.
    // aka, sum the dot product of each of the axes
    // https://en.wikipedia.org/wiki/Frobenius_inner_product
    // Another way to look at this is to just pairwise multiply every value in both arrays and sum the results.
    // since the rotations are "double sided", we need to negate one of the matrices and calculate the similarity again,
    // and return the larger of the two (the higher similarity value).
    // Another way we can do that is to just return the abs of the similarity.
    float fip = 0.0f;
    for (int i = 0; i < DIMENSION()*DIMENSION(); ++i)
        fip += A.rotation[i] * B.rotation[i];

    return std::fabsf(fip);
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_WhiteNoise(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    std::uniform_real_distribution<float> dist_angle(0.0f, 2.0f * c_pi);
    std::uniform_real_distribution<float> dist_offset(0.0f, 1.0f);

    for (PointHashData& p : hashDatas)
    {
        float angle = dist_angle(RNG());

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        p.angle = angle;

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = dist_offset(RNG());
    }
}

// -------------------------------------------------------------------------------

template <typename T, typename LAMBDA1, typename LAMBDA2>
void AlansGoodCandidateAlgorithm(std::vector<T>& results, int desiredItemCount, int candidateMultiplier, const LAMBDA2& GenerateRandomCandidate, const LAMBDA1& DifferenceScoreCalculator)
{
    // map candidate index to score
    struct CandidateScore
    {
        size_t index;
        float score;
    };
    typedef std::vector<CandidateScore> CandidateScores;
    static const size_t c_numScores = (1 << DIMENSION()) - 1;  // 2^(dimension)-1

    // make space for the results
    results.resize(desiredItemCount);

    // for each item we need to fill in
    for (int itemIndex = 0; itemIndex < desiredItemCount; ++itemIndex)
    {
        // calculate how many candidates we want to generate for this item
        int candidateCount = itemIndex * candidateMultiplier + 1;

        // generate the candidates
        std::vector<T> candidates;
        candidates.resize(candidateCount);
        for (T& candidate : candidates)
            candidate = GenerateRandomCandidate();

        // initialize the overall scores
        CandidateScores overallScores;
        overallScores.resize(candidateCount);
        for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
        {
            overallScores[candidateIndex].index = candidateIndex;
            overallScores[candidateIndex].score = 0.0f;
        }

        // allocate space for the individual scores
        CandidateScores scores;
        scores.resize(candidateCount);

        // score the candidates by each measure of scoring
        for (size_t scoreIndex = 0; scoreIndex < c_numScores; ++scoreIndex)
        {
            // for each candidate in this score index...
            for (size_t candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
            {
                const T& candidate = candidates[candidateIndex];

                // calculate the score of the candidate.
                // the score is the minimum distance to any other points
                float minimumDifferenceScore = FLT_MAX;
                for (int checkItemIndex = 0; checkItemIndex < itemIndex; ++checkItemIndex)
                {
                    float differenceScore = DifferenceScoreCalculator(candidate, results[checkItemIndex], scoreIndex);
                    minimumDifferenceScore = std::min(minimumDifferenceScore, differenceScore);
                }

                scores[candidateIndex].index = candidateIndex;
                scores[candidateIndex].score = minimumDifferenceScore;
            }

            // sort the scores from high to low
            std::sort(
                scores.begin(),
                scores.end(),
                [] (const CandidateScore& A, const CandidateScore& B)
                {
                    return A.score > B.score;
                }
            );

            // add the rank of this score a score for each candidate
            for (size_t candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
                overallScores[scores[candidateIndex].index].score += float(candidateIndex);
        }

        // sort the overall scores from low to high
        std::sort(
            overallScores.begin(),
            overallScores.end(),
            [] (const CandidateScore& A, const CandidateScore& B)
            {
                return A.score < B.score;
            }
        );

        // keep the point that had the lowest summed rank
        results[itemIndex] = candidates[overallScores[0].index];
    }
}

// -------------------------------------------------------------------------------

template <typename T, typename LAMBDA1, typename LAMBDA2>
void MitchelsBestCandidateAlgorithm (std::vector<T>& results, int desiredItemCount, int candidateMultiplier, const LAMBDA2& GenerateRandomCandidate, const LAMBDA1& DifferenceScoreCalculator)
{
    results.resize(desiredItemCount);

    // for each item we need to fill in
    for (int itemIndex = 0; itemIndex < desiredItemCount; ++itemIndex)
    {
        // calculate how many candidates we want to generate for this item
        int candidateCount = itemIndex * candidateMultiplier + 1;

        T bestCandidate;
        float bestCandidateMinimumDifferenceScore;

        // for each candidate
        for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
        {
            // make a randomized candidate
            T candidate = GenerateRandomCandidate();
            float minimumDifferenceScore = FLT_MAX;

            // the score of this candidate is the minimum difference from all existing items
            for (int checkItemIndex = 0; checkItemIndex < itemIndex; ++checkItemIndex)
            {
                float differenceScore = DifferenceScoreCalculator(candidate, results[checkItemIndex]);
                minimumDifferenceScore = std::min(minimumDifferenceScore, differenceScore);
            }

            // the candidate with the largest minimum distance is the one we want to keep
            if (candidateIndex == 0 || minimumDifferenceScore > bestCandidateMinimumDifferenceScore)
            {
                bestCandidate = candidate;
                bestCandidateMinimumDifferenceScore = minimumDifferenceScore;
            }
        }

        // keep the winning candidate
        results[itemIndex] = bestCandidate;
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_BlueNoise_1d(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::vector<float> angles, offsets;

    MitchelsBestCandidateAlgorithm(
        angles,
        HASHCOUNT(),
        5,
        [&]()
        {
            return dist(RNG());
        },
        [](const float& A, const float& B)
        {
            float diff = fabsf(B - A);
            if (diff > 0.5f)
                diff = 1.0f - diff;
            return diff;
        }
    );

    MitchelsBestCandidateAlgorithm(
        offsets,
        HASHCOUNT(),
        5,
        [&]()
        {
            return dist(RNG());
        },
        [](const float& A, const float& B)
        {
            float diff = fabsf(B - A);
            if (diff > 0.5f)
                diff = 1.0f - diff;
            return diff;
        }
    );

    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        float angle = angles[hashIndex] * c_pi;

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        hashDatas[hashIndex].angle = angle;

        hashDatas[hashIndex].rotation = { cosTheta, -sinTheta,
                                          sinTheta,  cosTheta };

        hashDatas[hashIndex].offsetX = offsets[hashIndex];
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_BlueNoise_2D(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::vector<TPoint> points;

    MitchelsBestCandidateAlgorithm(
        points,
        HASHCOUNT(),
        5,
        [&]()
        {
            TPoint ret;
            for (int i = 0; i < DIMENSION(); ++i)
                ret[i] = dist(RNG());
            return ret;
        },
        [](const TPoint& A, const TPoint& B)
        {
            float distSq = 0.0f;
            for (int i = 0; i < DIMENSION(); ++i)
            {
                float diff = fabsf(B[i] - A[i]);
                if (diff > 0.5f)
                    diff = 1.0f - diff;
                distSq += diff * diff;
            }
            return distSq;
        }
    );

    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        float angle = points[hashIndex][0] * c_pi;

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        hashDatas[hashIndex].angle = angle;

        hashDatas[hashIndex].rotation = { cosTheta, -sinTheta,
                                          sinTheta,  cosTheta };

        hashDatas[hashIndex].offsetX = points[hashIndex][1];
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_BlueNoise_2DProjective(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::vector<TPoint> points;

    static const size_t c_numScores = (1 << DIMENSION()) - 1;
    TPoint axisMasks[c_numScores];
    for (size_t i = 1; i <= c_numScores; ++i)
    {
        TPoint& axisMask = axisMasks[i - 1];

        for (size_t j = 0; j < DIMENSION(); ++j)
        {
            if (i & (size_t(1) << j))
                axisMask[j] = 1.0f;
            else
                axisMask[j] = 0.0f;
        }
    }

    AlansGoodCandidateAlgorithm(
        points,
        HASHCOUNT(),
        20,
        [&]()
        {
            TPoint ret;
            for (int i = 0; i < DIMENSION(); ++i)
                ret[i] = dist(RNG());
            return ret;
        },
        [&](const TPoint& A, const TPoint& B, size_t scoreIndex)
        {
            const TPoint& axisMask = axisMasks[scoreIndex];

            float distSq = 0.0f;
            for (int i = 0; i < DIMENSION(); ++i)
            {
                float diff = fabsf(B[i] - A[i]) * axisMask[i];
                if (diff > 0.5f)
                    diff = 1.0f - diff;
                distSq += diff * diff;
            }
            return sqrtf(distSq);
        }
    );

    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        float angle = points[hashIndex][0] * c_pi;

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        hashDatas[hashIndex].angle = angle;

        hashDatas[hashIndex].rotation = { cosTheta, -sinTheta,
                                          sinTheta,  cosTheta };

        hashDatas[hashIndex].offsetX = points[hashIndex][1];
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_Uniform(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    int columns = int(ceilf(sqrtf(float(HASHCOUNT()))));
    int rows = int(ceilf(float(HASHCOUNT()) / float(columns)));

    for (int i = 0; i < HASHCOUNT(); ++i)
    {
        PointHashData& p = hashDatas[i];

        int x = i % columns;
        int y = i / columns;

        float percentX = (float(x) + 0.5f) / float(columns);
        float percentY = (float(y) + 0.5f) / float(rows);

        // choose an angle from 0 to 180 degrees, because angle > 180 degrees is redundant.
        float angle = percentX * c_pi;

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        p.angle = angle;

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = percentY;
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_GoldenRatio(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    float initialValueAngle = dist(RNG());
    float initialValueOffset = dist(RNG());

    for (int i = 0; i < HASHCOUNT(); ++i)
    {
        PointHashData& p = hashDatas[i];

        // choose an angle from 0 to 180 degrees, because angle > 180 degrees is redundant and actually harmful.
        // the angles are "double sided".
        // check the angle test images for more information!
        float angle = std::fmodf(initialValueAngle + float(i) * c_goldenRatioConjugate, 1.0f) * c_pi;

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        p.angle = angle;

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = std::fmodf(initialValueOffset + float(i) * c_goldenRatioConjugate, 1.0f);
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_GoldenRatioGeneralized(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    float initialValue = dist(RNG());

    for (int i = 0; i < HASHCOUNT(); ++i)
    {
        PointHashData& p = hashDatas[i];

        // choose an angle from 0 to 180 degrees, because angle > 180 degrees is redundant and actually harmful.
        // the angles are "double sided".
        // check the angle test images for more information!
        float angle = std::fmodf(initialValue + float(i) * c_goldenRatio2Conjugate, 1.0f) * c_pi;

        float cosTheta = std::cosf(angle);
        float sinTheta = std::sinf(angle);

        p.angle = angle;

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = std::fmodf(initialValue + float(i) * c_goldenRatio2Conjugate * c_goldenRatio2Conjugate, 1.0f);
    }
}

// -------------------------------------------------------------------------------
void SaveImage(const char* fileName, Image& image)
{
    image.m_pixels[((image.m_width*image.m_height - 1) * 4) + 3] = 0; // make the last pixel be transparent so eg twitter doesn't use jpg compression.
    stbi_write_png(fileName, image.m_width, image.m_height, 4, image.m_pixels.data(), 0);
}

// -------------------------------------------------------------------------------

float SmoothStep(float value, float min, float max)
{
    float x = (value - min) / (max - min);
    x = std::min(x, 1.0f);
    x = std::max(x, 0.0f);

    return 3.0f * x * x - 2.0f * x * x * x;
}

// -------------------------------------------------------------------------------

template <typename T>
T Lerp(T A, T B, float t)
{
    return T(float(A) * (1.0f - t) + float(B) * t);
}

// -------------------------------------------------------------------------------

void DrawLine(Image& image, int x1, int y1, int x2, int y2, uint8 R, uint8 G, uint8 B)
{
    // pad the AABB of pixels we scan, to account for anti aliasing
    int startX = std::max(std::min(x1, x2) - 4, 0);
    int startY = std::max(std::min(y1, y2) - 4, 0);
    int endX = std::min(std::max(x1, x2) + 4, image.m_width - 1);
    int endY = std::min(std::max(y1, y2) + 4, image.m_height - 1);

    // if (x1,y1) is A and (x2,y2) is B, get a normalized vector from A to B called AB
    float ABX = float(x2 - x1);
    float ABY = float(y2 - y1);
    float ABLen = std::sqrtf(ABX*ABX + ABY * ABY);
    ABX /= ABLen;
    ABY /= ABLen;

    // scan the AABB of our line segment, drawing pixels for the line, as is appropriate
    for (int iy = startY; iy <= endY; ++iy)
    {
        uint8* pixel = &image.m_pixels[(iy * image.m_width + startX) * 4];
        for (int ix = startX; ix <= endX; ++ix)
        {
            // project this current pixel onto the line segment to get the closest point on the line segment to the point
            float ACX = float(ix - x1);
            float ACY = float(iy - y1);
            float lineSegmentT = ACX * ABX + ACY * ABY;
            lineSegmentT = std::min(lineSegmentT, ABLen);
            lineSegmentT = std::max(lineSegmentT, 0.0f);
            float closestX = float(x1) + lineSegmentT * ABX;
            float closestY = float(y1) + lineSegmentT * ABY;

            // calculate the distance from this pixel to the closest point on the line segment
            float distanceX = float(ix) - closestX;
            float distanceY = float(iy) - closestY;
            float distance = std::sqrtf(distanceX*distanceX + distanceY * distanceY);

            // use the distance to figure out how transparent the pixel should be, and apply the color to the pixel
            float alpha = SmoothStep(distance, 2.0f, 0.0f);

            if (alpha > 0.0f)
            {
                pixel[0] = Lerp(pixel[0], R, alpha);
                pixel[1] = Lerp(pixel[1], G, alpha);
                pixel[2] = Lerp(pixel[2], B, alpha);
            }

            pixel += 4;
        }
    }
}

// -------------------------------------------------------------------------------

void ClearImage(Image& image, uint8 R, uint8 G, uint8 B)
{
    uint8* pixel = image.m_pixels.data();
    for (int i = 0, c = image.m_width * image.m_height; i < c; ++i)
    {
        pixel[0] = R;
        pixel[1] = G;
        pixel[2] = B;
        pixel[3] = 255;
        pixel += 4;
    }
}

// -------------------------------------------------------------------------------

void DrawCircle(Image& image, int cx, int cy, int radius, uint8 R, uint8 G, uint8 B)
{
    int startX = std::max(cx - radius - 4, 0);
    int startY = std::max(cy - radius - 4, 0);
    int endX = std::min(cx + radius + 4, image.m_width - 1);
    int endY = std::min(cy + radius + 4, image.m_height - 1);

    for (int iy = startY; iy <= endY; ++iy)
    {
        float dy = float(cy - iy);
        uint8* pixel = &image.m_pixels[(iy * image.m_width + startX) * 4];
        for (int ix = startX; ix <= endX; ++ix)
        {
            float dx = float(cx - ix);

            float distance = std::max(std::sqrtf(dx * dx + dy * dy) - float(radius), 0.0f);

            float alpha = SmoothStep(distance, 2.0f, 0.0f);

            if (alpha > 0.0f)
            {
                pixel[0] = Lerp(pixel[0], R, alpha);
                pixel[1] = Lerp(pixel[1], G, alpha);
                pixel[2] = Lerp(pixel[2], B, alpha);
            }

            pixel += 4;
        }
    }
}

// -------------------------------------------------------------------------------

void DrawHashBuckets(Image& image, float cellSize, const PointHashData& pointHashData, uint8 R, uint8 G, uint8 B)
{
    uint8* pixel = image.m_pixels.data();
    for (int iy = 0; iy < image.m_height; ++iy)
    {
        float posY = Lerp(-float(POINTDOMAIN()), float(POINTDOMAIN()), float(iy) / float(image.m_height));
        for (int ix = 0; ix < image.m_width; ++ix)
        {
            float posX = Lerp(-float(POINTDOMAIN()), float(POINTDOMAIN()), float(ix) / float(image.m_width));

            // transform pixel into the hash bucket space
            TPoint rawPoint = {float(posX), float(posY)};
            TPoint rotatedPoint = Multiply(rawPoint, pointHashData.rotation);
            rotatedPoint[0] += pointHashData.offsetX;

            float distX = std::fmodf(fabsf(rotatedPoint[0]), 1.0f);
            float distance = 0.5f - fabsf(distX - 0.5f);

            // draw the hash bucket lines
            float alpha = SmoothStep(distance, 4.0f * float(POINTDOMAIN()) / float(IMAGESIZE()), 0.0f);

            if (alpha > 0.0f)
            {
                pixel[0] = Lerp(pixel[0], R, alpha);
                pixel[1] = Lerp(pixel[1], G, alpha);
                pixel[2] = Lerp(pixel[2], B, alpha);
            }
            pixel += 4;
        }
    }
}

// -------------------------------------------------------------------------------

void ReportQueryAsImage(const LHS& lhs, const TPoint& queryPoint, const std::unordered_map<PointID, int>& results, const std::vector<TPoint>& points, const char* name)
{
    char fileName[256];
    sprintf_s(fileName, "out/%s.png", name);

    int circleRadius = int(float(IMAGESIZE()) / 100.0f);
    int colorCircleRadius = std::min(int(float(IMAGESIZE()) / 100.0f * 0.9f), circleRadius - 1);

    // create the image and fill it with white.
    Image imageTop(IMAGESIZE(), IMAGESIZE());

    // draw the hash buckets
    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        const PointHashData& pointHashData = lhs.GetPointHashData(hashIndex);
        DrawHashBuckets(imageTop, float(IMAGESIZE()) * 0.5f / float(POINTDOMAIN()), pointHashData, 192, 192, 192);
    }

    // draw all the points as black dots
    for (TPoint p : points)
    {
        for (float& f : p)
        {
            f /= (2.0f * float(POINTDOMAIN()));
            f += 0.5f;
        }

        int x = int(p[0] * float(IMAGESIZE() - 1) + 0.5f);
        int y = int(p[1] * float(IMAGESIZE() - 1) + 0.5f);

        DrawCircle(imageTop, x, y, circleRadius, 0, 0, 0);
    }

    // draw the results based on their match count
    for (auto it : results)
    {
        if (it.second < MINHASHCOUNT())
            continue;

        TPoint p = points[(uint32)it.first];

        for (float& f : p)
        {
            f /= (2.0f * float(POINTDOMAIN()));
            f += 0.5f;
        }

        int x = int(p[0] * float(IMAGESIZE() - 1) + 0.5f);
        int y = int(p[1] * float(IMAGESIZE() - 1) + 0.5f);

        float percentMatched = (float(it.second) - float(MINHASHCOUNT())) / float(HASHCOUNT() - MINHASHCOUNT());
        uint8 color = uint8(percentMatched * 255.0f);

        DrawCircle(imageTop, x, y, colorCircleRadius, 255, color, 0);
    }

    // draw the query point
    {
        TPoint p = queryPoint;
        for (float& f : p)
        {
            f /= (2.0f * float(POINTDOMAIN()));
            f += 0.5f;
        }
        int x = int(p[0] * float(IMAGESIZE() - 1) + 0.5f);
        int y = int(p[1] * float(IMAGESIZE() - 1) + 0.5f);

        DrawCircle(imageTop, x, y, circleRadius, 0, 0, 255);
    }

    // draw the footer - show the angle and offset distribution on a number line
    Image imageBottom(IMAGESIZE(), IMAGESIZE());
    Image imageSamples(50, 50);  // TODO: make this be a #define at the top for how large it is
    ImageComplex imageSamplesDFTComplex(imageSamples.m_width, imageSamples.m_height);
    DrawLine(imageBottom, 0, 0, IMAGESIZE(), 0, 0, 0, 0);
    DrawLine(imageBottom, 0, IMAGESIZE()-1, IMAGESIZE(), IMAGESIZE()-1, 0, 0, 0);

    int graphSize = int(float(IMAGESIZE()) * 0.7f);
    int padding = (IMAGESIZE() - graphSize) / 2;

    DrawLine(imageBottom, padding, padding, padding + graphSize, padding, 128, 128, 255);
    DrawLine(imageBottom, padding, padding + graphSize, padding + graphSize, padding + graphSize, 128, 128, 255);

    DrawLine(imageBottom, padding, padding, padding, padding + graphSize, 128, 128, 255);
    DrawLine(imageBottom, padding + graphSize, padding, padding + graphSize, padding + graphSize, 128, 128, 255);

    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        const PointHashData& p = lhs.GetPointHashData(hashIndex);
        float hashPercent = float(hashIndex) / float(HASHCOUNT() - 1);
        uint8 color = uint8(255.0f * hashPercent);

        float anglePercent = std::fmodf(p.angle / c_pi, 1.0f);
        int dotX = padding + int(anglePercent * float(graphSize));

        int dotY = padding + int(p.offsetX * float(graphSize));

        // draw the 2d dots
        DrawCircle(imageBottom, dotX, dotY, 2, 0, color, 0);

        // draw the 1d angle dots to the right
        DrawCircle(imageBottom, padding + graphSize + padding / 2, dotY, 2, 0, color, 0);

        // draw the 1d offset dots below
        DrawCircle(imageBottom, dotX, padding + graphSize + padding / 2, 2, 0, color, 0);

        // make the sample image
        {
            int x = int(anglePercent * float(imageSamples.m_width - 1) + 0.5f);
            int y = int(p.offsetX  * float(imageSamples.m_height - 1) + 0.5f);
            uint8* imageSamplePixel = &imageSamples.m_pixels[(y*imageSamples.m_width + x) * 4];
            imageSamplePixel[0] = 0;
            imageSamplePixel[1] = 0;
            imageSamplePixel[2] = 0;
            imageSamplePixel[3] = 255;
        }
    }

    // DFT the image samples
    DFTImage(imageSamples, imageSamplesDFTComplex);

    // combine the images
    Image image(IMAGESIZE(), IMAGESIZE() * 2);
    memcpy(image.m_pixels.data(), imageTop.m_pixels.data(), imageTop.m_pixels.size());
    memcpy(&image.m_pixels[IMAGESIZE()*IMAGESIZE() * 4], imageBottom.m_pixels.data(), imageBottom.m_pixels.size());

    // TODO: formalize an "append vertical" function and copy row by row into the main image somewhere.
    // TODO: make the dft functions (getmag etc) go y on outside loop, x on inside for cache
    Image imageSamplesDFT(imageSamples.m_width, imageSamples.m_height);
    GetMagnitudeData(imageSamplesDFTComplex, imageSamplesDFT);
    //SaveImage(fileName, imageSamplesDFT);

    SaveImage(fileName, image);
}

// -------------------------------------------------------------------------------

void ReportQuery (const LHS& lhs, const TPoint& queryPoint, const std::vector<TPoint>& points, const char* name)
{
    printf("\n=====%s=====\n", name);

    std::unordered_map<PointID, int> results;
    lhs.Query(queryPoint, results);

    for (int hashMatchCount = HASHCOUNT(); hashMatchCount >= MINHASHCOUNT(); --hashMatchCount)
    {
        bool foundMatch = false;
        for (auto it = results.begin(); it != results.end(); ++it)
        {
            if (it->second == hashMatchCount)
            {
                if (!foundMatch)
                {
                    printf("\nHash Match %i\n", hashMatchCount);
                    foundMatch = true;
                }

                const TPoint& point = lhs.GetPoint(it->first);
                printfPoint(nullptr, point);
            }
        }
    }

    ReportQueryAsImage(lhs, queryPoint, results, points, name);
}

// -------------------------------------------------------------------------------

void AngleImageTest()
{
    static const int c_numAngles = 10;

    Image imageAngle(IMAGESIZE(), IMAGESIZE());
    DrawCircle(imageAngle, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(imageAngle, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    Image imageAngleDoubleSided(IMAGESIZE(), IMAGESIZE());
    DrawCircle(imageAngleDoubleSided, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(imageAngleDoubleSided, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    Image imageHalfAngle(IMAGESIZE(), IMAGESIZE());
    DrawCircle(imageHalfAngle, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(imageHalfAngle, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    Image imageHalfAngleDoubleSided(IMAGESIZE(), IMAGESIZE());
    DrawCircle(imageHalfAngleDoubleSided, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(imageHalfAngleDoubleSided, IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    for (int i = 0; i < c_numAngles; ++i)
    {
        // use golden ratio to make angles between 0 and 2*pi
        {
            float angleRad = std::fmodf(float(i) * c_goldenRatioConjugate, 1.0f) * 2.0f * c_pi;
            float angleDeg = angleRad * 180.0f / c_pi;

            int lineEndX = int(cos(angleRad) * (float(IMAGESIZE() / 2) - 4) + float(IMAGESIZE() / 2));
            int lineEndY = int(sin(angleRad) * (float(IMAGESIZE() / 2) - 4) + float(IMAGESIZE() / 2));

            lineEndY = IMAGESIZE() - lineEndY;

            uint8 color = uint8(255.0f * float(i) / float(c_numAngles));

            DrawLine(imageAngle, IMAGESIZE() / 2, IMAGESIZE() / 2, lineEndX, lineEndY, color, color, color);

            // also draw the double sided version
            DrawLine(imageAngleDoubleSided, IMAGESIZE() - lineEndX, IMAGESIZE() - lineEndY, lineEndX, lineEndY, color, color, color);
        }

        // use golden ratio to make angles between 0 and pi
        {
            float angleRad = std::fmodf(float(i) * c_goldenRatioConjugate, 1.0f) * c_pi;
            float angleDeg = angleRad * 180.0f / c_pi;

            int lineEndX = int(cos(angleRad) * (float(IMAGESIZE() / 2) - 4) + float(IMAGESIZE() / 2));
            int lineEndY = int(sin(angleRad) * (float(IMAGESIZE() / 2) - 4) + float(IMAGESIZE() / 2));

            lineEndY = IMAGESIZE() - lineEndY;

            uint8 color = uint8(255.0f * float(i) / float(c_numAngles));

            DrawLine(imageHalfAngle, IMAGESIZE() / 2, IMAGESIZE() / 2, lineEndX, lineEndY, color, color, color);

            // also draw the double sided version
            DrawLine(imageHalfAngleDoubleSided, IMAGESIZE() - lineEndX, IMAGESIZE() - lineEndY, lineEndX, lineEndY, color, color, color);
        }
    }

    SaveImage("out/AngleTest.png", imageAngle);
    SaveImage("out/AngleTestDS.png", imageAngleDoubleSided);
    SaveImage("out/AngleTestHalf.png", imageHalfAngle);
    SaveImage("out/AngleTestHalfDS.png", imageHalfAngleDoubleSided);
}

// -------------------------------------------------------------------------------

int main(int argc, char** argv)
{
#if DOANGLEIMAGETEST()
    AngleImageTest();
#endif

    std::uniform_real_distribution<float> dist(-float(POINTDOMAIN()), float(POINTDOMAIN()));

    // make the random point set
    std::vector<TPoint> points;
    points.resize(NUMPOINTS());
    for (TPoint& p : points)
    {
        for (float& f : p)
            f = dist(RNG());
    }

    // add the points to the locality sensitive hashing
    LHS lhs_whiteNoise(points, GeneratePointHashDatas_WhiteNoise);
    LHS lhs_blueNoise_1d(points, GeneratePointHashDatas_BlueNoise_1d);
    LHS lhs_blueNoise_2D(points, GeneratePointHashDatas_BlueNoise_2D);
    LHS lhs_blueNoise_2DProjective(points, GeneratePointHashDatas_BlueNoise_2DProjective);
    LHS lhs_uniform(points, GeneratePointHashDatas_Uniform);
    LHS lhs_goldenRatio(points, GeneratePointHashDatas_GoldenRatio);
    LHS lhs_goldenRatioGeneralized(points, GeneratePointHashDatas_GoldenRatioGeneralized);

    // do some queries
    for (int i = 0; i < 1; ++i)
    {
        TPoint queryPoint;
        for (float& f : queryPoint)
            f = dist(RNG());

        printfPoint("Query Point", queryPoint);

        ReportQuery(lhs_whiteNoise, queryPoint, points, "whiteNoise");
        ReportQuery(lhs_blueNoise_1d, queryPoint, points, "blueNoise_1d");
        ReportQuery(lhs_blueNoise_2D, queryPoint, points, "blueNoise_2d");
        ReportQuery(lhs_blueNoise_2DProjective, queryPoint, points, "blueNoise_2dProjective");
        ReportQuery(lhs_uniform, queryPoint, points, "uniform");
        ReportQuery(lhs_goldenRatio, queryPoint, points, "goldenRatio");
        ReportQuery(lhs_goldenRatioGeneralized, queryPoint, points, "goldenRatioGeneralized");
    }

    // TODO: ground truth, how? maybe visually show points missed or something? could calculate std deviations as concentric rings.

    return 0;
}

/*
 TODO:
 * what is the gap in the angular wrap around? is it 2x too large?
 * DFT of 2d dots and 1d projections.
? are the buckets too small? i think they might be... instead of having the domain thing for points, could have a bucket size scale.
* maybe an option to color cells based on how many hash collisions they have with the query point
* tyler's article says there is an optimal for white noise. check out what that is, and try it? maybe compare things with that optimal value?
* there is likely a point where white noise beats out blue / LDS. find it?
* maybe try halton or something
* maybe try that "low discrepancy blue noise" thing?
* check through the code. maybe delete things you don't need naymore?
* find out how to do n Owen-scrambled Sobol. maybe would be the ideal thing?
* random matrices if needed for higher dimensions, from tyler: One method is to choose normal randoms per matrix entry, take the QR decomp of that, and now Q is a random rotation matrix, essentially.

Tests:
* white noise
* uniform
* blue noise
* LDS (golden ratio?) - yes. generalized golden ratio. http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
* compare vs linear search to get ground truth

* Some way to gather stats, but also want to make diagrams i think. Maybe some animated like the other blog post does?
* and maybe calculate mean and variance of cell size if you can

* higher dimensional testing? to see how blue noise / golden ratio / etc behave



----- LATER -----

* maybe a second blog post on bloom filter?
 * blue noise is about quality of results.
 * bloom filter is about performance, but may also affect quality
 * the paper from tyler sort of mentions a bloom filter related thing: "if we�re willing to probabilistically verify hash keys at the cost of some false collisions, we can still achieve O(d log d) lookups" 
 * basically, if you get hashed items, you need to resolve them to points which takes longer. if you work only w/ hashes you can get false positives but it's faster. A bit bloom filter ish.

* maybe do 1d integration and compare blue noise, white noise, lds?
 * to show that blue noise has same error decay as white noise but starts lower?
 * could compare to owen scrambled sobol, and low discrepancy blue noise sampling.

Links:
- inspiration: http://unboxresearch.com/articles/lsh_post1.html   aka http://unboxresearch.com/articles/lsh_post1.html
* another way to do it: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.215.7690&rep=rep1&type=pdf
 * make regular cells of any shape. store a list at each cell corner which has each point that is in a cell that uses that corner.
 * talks about simplices fighting the curse of dimensionality.
 * regular simplices can't tile higher dimensions though so it has an alternate formulation.
- maybe some blue noise articles?
* frobenius inner product for checking similarity of rotation matrices: https://en.wikipedia.org/wiki/Frobenius_inner_product
* random rotation matrices: https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices
* lattice based LHSs are apparently optimal: https://arxiv.org/pdf/1712.08558.pdf
 * but they can totally be randomized still! do a random rotation / offset on the point before feeding it into the deterministic lattice.
* projective blue noise: http://resources.mpi-inf.mpg.de/ProjectiveBlueNoise/ProjectiveBlueNoise.pdf

Notes:
* have a white noise [0,1] random value to start golden ratio by. You get the benefits of random values i believe, like you do with white noise.
 * do you get all of them?
* Does it make sense to do golden ratio between 0 and pi instead of 0 and 2 pi? yes. check the image tests. The angles are "double sided".
* mitchell's best candidate algorithm for generating blue noise has a tuneable parameter. Too low and you get white noise. Too high and you get regular sampling. not super great for getting good results out of the box.
* generalized golden ratio doesn't look that LDS at lower hash counts, but it does at higher hash counts.  Maybe not useful since hash count should (?) always be same as # of faces on a simplex for that dimension?
* projective blue noise: I'm counting single axis distances to be worth as much as 1/sqrt(2) the 2d distance.
 * if i were doing 3d points, and looking at a 2d subspace i would count the 2d distances as sqrt(2)/sqrt(3).
 * if i were doing 3d points, and looking at a 1d subspace, i would count the 1d distances as 1/sqrt(3).
 * the projective blue noise paper uses sphere packing for specific number of points, which doesn't quite translate to best candidate.
* My "good score algorithm" grows in complexity very quickly for the number of dimensions


----- Good Candidate algorithm -----

co-author with brandon!

Motivation:
  Dart throwing and best candidate make progressive sequences, but best candidate gives better results.
  Actually, dart throwing isn't progressive out of the box really! You have a radius!
  Anyways, not as high quality as voronoi relaxation method, but it's progressive which is neat. Maybe eg 1024 points as constants in a shader, use however many you want, subset of that.

Do the tests from the multijittered sampling: disk, step, gaussian, bilinear.

compare to... ??
* whtie noise
* blue noise
* owen scrambled sobol
? projective blue noise sampling via other methods?
? low discrepancy blue noise?

other links:
renderman: https://graphics.pixar.com/library/RendermanTog2018/paper.pdf
projective blue noise article: http://resources.mpi-inf.mpg.de/ProjectiveBlueNoise/
multijittered sampling: https://graphics.pixar.com/library/ProgressiveMultiJitteredSampling/

----- LANDFILL -----

* do we need the full matrix? i don't think so, if we are only takiung the resulting x component!
 * maybe we don't care. implementation detail (optimization)

*/