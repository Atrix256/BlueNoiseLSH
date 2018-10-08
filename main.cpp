// Configurable settings
#define DIMENSION()    2      // dimensionality of the data points
#define NUMPOINTS()    1000   // the number of randomly generated (white noise) data points
#define HASHCOUNT()    8      // how many hashes are done for each point
#define MINHASHCOUNT() 2      // the minimum number of hash collisions that must occur for a point to be considered at all close
#define POINTDOMAIN()  5      // the coordinates go from - this value to + this value

#define IMAGESIZE()       500 // the size of the image - width and height both.  T
#define IMAGEFOOTERSIZE() 500 // the size of the "footer" under the image that shows distribution of angle and offset

#define DOANGLEIMAGETEST() 1

// TODO: tune HASHCOUNT() and the other parameters!

#include <random>
#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <stdint.h>

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
    static_assert(DIMENSION() == 2, "This function only works with 2d rotation matrices");

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

void GeneratePointHashDatas_BlueNoise(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    static_assert(DIMENSION() == 2, "This function only works with 2d rotation matrices");

    std::uniform_real_distribution<float> dist_angle(0.0f, 2.0f * c_pi);
    std::uniform_real_distribution<float> dist_offset(0.0f, 1.0f);

    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        int numCandidates = hashIndex * 10 + 1;

        PointHashData bestPointHashData;
        float bestPointHashDataScore;

        for (int candidateIndex = 0; candidateIndex < numCandidates; ++candidateIndex)
        {
            PointHashData candidatePointHashData;

            float angle = dist_angle(RNG());

            float cosTheta = std::cosf(angle);
            float sinTheta = std::sinf(angle);

            candidatePointHashData.angle = angle;

            candidatePointHashData.rotation = { cosTheta, -sinTheta,
                                                sinTheta,  cosTheta };

            // TODO: we need to find a way to weigh offset in the mix for SimilarityScore(). after that, make this work.
            candidatePointHashData.offsetX = dist_offset(RNG());

            // the score of this candidate is the most similar it is to any existing data point
            float maxScore = 0.0f;
            for (int dataIndex = 0; dataIndex < hashIndex; ++dataIndex)
            {
                float score = SimilarityScore(candidatePointHashData, hashDatas[dataIndex]);
                if (dataIndex == 0)
                    maxScore = score;
                else
                    maxScore = std::max(score, maxScore);
            }

            // the best candidate will have the smallest score, meaning it's least similar to all existing data points
            if (candidateIndex == 0 || maxScore < bestPointHashDataScore)
            {
                bestPointHashData = candidatePointHashData;
                bestPointHashDataScore = maxScore;
            }
        }

        // keep the best candidate found
        hashDatas[hashIndex] = bestPointHashData;
    }
}

// -------------------------------------------------------------------------------

void GeneratePointHashDatas_Uniform(std::array<PointHashData, HASHCOUNT()>& hashDatas)
{
    static_assert(DIMENSION() == 2, "This function only works with 2d rotation matrices");

    // To sample evenly across the two dimensional space of angle and offset, we are going to treat it as a square
    // of dimensions (1,1) that we need to place HASHCOUNT() points in uniformly.
    // So, we need to calculate the number of rows and columns we'll need.
    // The number after each square (power of 2) number of points is when a new row is going to be added.
    // So, rows == columns == ceil(sqrt(HASHCOUNT()).
    // This would generalize to higher dimensions, replacing sqrt with the appropriate root

    // TODO: how to deal with the remainder? I guess you could choose a row or column and re-center them on that axis. is that most uniform?
    // TODO: do we want to treat each axis equally? if not, would need to adjust the calculation somehow.

    int size = int(ceilf(sqrtf(float(HASHCOUNT()))));

    for (int i = 0; i < HASHCOUNT(); ++i)
    {
        PointHashData& p = hashDatas[i];

        int x = i % size;
        int y = i / size;

        float percentX = (float(x) + 0.5f) / float(size);
        float percentY = (float(y) + 0.5f) / float(size);

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
    static_assert(DIMENSION() == 2, "This function only works with 2d rotation matrices");

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
    static_assert(DIMENSION() == 2, "This function only works with 2d rotation matrices");

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
void SaveImage(const char* fileName, int width, int height, uint8* pixels)
{
    pixels[((width*height - 1) * 4) + 3] = 0; // make the last pixel be transparent so eg twitter doesn't use jpg compression.
    stbi_write_png(fileName, width, height, 4, pixels, 0);
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

void DrawLine(std::vector<uint8>& pixels, int imageWidth, int imageHeight, int x1, int y1, int x2, int y2, uint8 R, uint8 G, uint8 B)
{
    // pad the AABB of pixels we scan, to account for anti aliasing
    int startX = std::max(std::min(x1, x2) - 4, 0);
    int startY = std::max(std::min(y1, y2) - 4, 0);
    int endX = std::min(std::max(x1, x2) + 4, imageWidth - 1);
    int endY = std::min(std::max(y1, y2) + 4, imageHeight - 1);

    // if (x1,y1) is A and (x2,y2) is B, get a normalized vector from A to B called AB
    float ABX = float(x2 - x1);
    float ABY = float(y2 - y1);
    float ABLen = std::sqrtf(ABX*ABX + ABY * ABY);
    ABX /= ABLen;
    ABY /= ABLen;

    // scan the AABB of our line segment, drawing pixels for the line, as is appropriate
    for (int iy = startY; iy <= endY; ++iy)
    {
        uint8* pixel = &pixels[(iy * imageWidth + startX) * 4];
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

void ClearImage(std::vector<uint8>& pixels, int imageWidth, int imageHeight, uint8 R, uint8 G, uint8 B)
{
    uint8* pixel = pixels.data();
    for (int i = 0, c = imageWidth * imageHeight; i < c; ++i)
    {
        pixel[0] = R;
        pixel[1] = G;
        pixel[2] = B;
        pixel[3] = 255;
        pixel += 4;
    }
}

// -------------------------------------------------------------------------------

void DrawCircle(std::vector<uint8>& pixels, int imageWidth, int imageHeight, int cx, int cy, int radius, uint8 R, uint8 G, uint8 B)
{
    int startX = std::max(cx - radius - 4, 0);
    int startY = std::max(cy - radius - 4, 0);
    int endX = std::min(cx + radius + 4, imageWidth - 1);
    int endY = std::min(cy + radius + 4, imageHeight - 1);

    for (int iy = startY; iy <= endY; ++iy)
    {
        float dy = float(cy - iy);
        uint8* pixel = &pixels[(iy * imageWidth + startX) * 4];
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

void DrawHashBuckets(std::vector<uint8>& pixels, int imageWidth, int imageHeight, float cellSize, const PointHashData& pointHashData, uint8 R, uint8 G, uint8 B)
{
    uint8* pixel = pixels.data();
    for (int iy = 0; iy < imageHeight; ++iy)
    {
        float posY = Lerp(-float(POINTDOMAIN()), float(POINTDOMAIN()), float(iy) / float(imageHeight));
        for (int ix = 0; ix < imageWidth; ++ix)
        {
            float posX = Lerp(-float(POINTDOMAIN()), float(POINTDOMAIN()), float(ix) / float(imageHeight));

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
    static_assert(DIMENSION() == 2, "This function only works with 2d points");

    char fileName[256];
    sprintf_s(fileName, "out/%s.png", name);

    int circleRadius = int(float(IMAGESIZE()) / 100.0f);
    int colorCircleRadius = std::min(int(float(IMAGESIZE()) / 100.0f * 0.9f), circleRadius - 1);

    // create the image and fill it with white.
    std::vector<uint8> pixels;
    pixels.resize(IMAGESIZE()*(IMAGESIZE()+IMAGEFOOTERSIZE()) * 4);
    std::fill(pixels.begin(), pixels.end(), 255);

    // draw the hash buckets
    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        const PointHashData& pointHashData = lhs.GetPointHashData(hashIndex);
        DrawHashBuckets(pixels, IMAGESIZE(), IMAGESIZE(), float(IMAGESIZE()) * 0.5f / float(POINTDOMAIN()), pointHashData, 192, 192, 192);
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

        DrawCircle(pixels, IMAGESIZE(), IMAGESIZE(), x, y, circleRadius, 0, 0, 0);
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

        DrawCircle(pixels, IMAGESIZE(), IMAGESIZE(), x, y, colorCircleRadius, 255, color, 0);
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

        DrawCircle(pixels, IMAGESIZE(), IMAGESIZE(), x, y, circleRadius, 0, 0, 255);
    }

    // draw the footer - show the angle and offset distribution on a number line
    DrawLine(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), 0, IMAGESIZE(), IMAGESIZE(), IMAGESIZE(), 0, 0, 0);

    int graphSize = int(float(IMAGEFOOTERSIZE()) * 0.7f);
    int padding = (IMAGEFOOTERSIZE() - graphSize) / 2;

    DrawLine(pixels, IMAGESIZE(), IMAGESIZE()+IMAGEFOOTERSIZE(), padding, IMAGESIZE() + padding, padding + graphSize, IMAGESIZE() + padding, 128, 128, 255);
    DrawLine(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), padding, IMAGESIZE() + padding + graphSize, padding + graphSize, IMAGESIZE() + padding + graphSize, 128, 128, 255);

    DrawLine(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), padding, IMAGESIZE() + padding, padding, IMAGESIZE() + padding + graphSize, 128, 128, 255);
    DrawLine(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), padding + graphSize, IMAGESIZE() + padding, padding + graphSize, IMAGESIZE() + padding + graphSize, 128, 128, 255);

    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        const PointHashData& p = lhs.GetPointHashData(hashIndex);
        float hashPercent = float(hashIndex) / float(HASHCOUNT() - 1);
        uint8 color = uint8(255.0f * hashPercent);

        float anglePercent = std::fmodf(p.angle / c_pi, 1.0f);
        int dotX = padding + int(anglePercent * float(graphSize));

        int dotY = IMAGESIZE() + padding + int(p.offsetX * float(graphSize));

        // draw the 2d dots
        DrawCircle(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), dotX, dotY, 2, 0, color, 0);

        // draw the 1d angle dots to the right
        DrawCircle(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), padding + graphSize + padding / 2, dotY, 2, 0, color, 0);

        // draw the 1d offset dots below
        DrawCircle(pixels, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), dotX, IMAGESIZE() + padding + graphSize + padding / 2, 2, 0, color, 0);
    }

    SaveImage(fileName, IMAGESIZE(), IMAGESIZE() + IMAGEFOOTERSIZE(), pixels.data());
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

    std::vector<uint8> pixelsAngle;
    pixelsAngle.resize(IMAGESIZE()*IMAGESIZE() * 4);
    std::fill(pixelsAngle.begin(), pixelsAngle.end(), 255);
    DrawCircle(pixelsAngle, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(pixelsAngle, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    std::vector<uint8> pixelsAngleDoubleSided;
    pixelsAngleDoubleSided.resize(IMAGESIZE()*IMAGESIZE() * 4);
    std::fill(pixelsAngleDoubleSided.begin(), pixelsAngleDoubleSided.end(), 255);
    DrawCircle(pixelsAngleDoubleSided, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(pixelsAngleDoubleSided, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    std::vector<uint8> pixelsHalfAngle;
    pixelsHalfAngle.resize(IMAGESIZE()*IMAGESIZE() * 4);
    std::fill(pixelsHalfAngle.begin(), pixelsHalfAngle.end(), 255);
    DrawCircle(pixelsHalfAngle, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(pixelsHalfAngle, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

    std::vector<uint8> pixelsHalfAngleDoubleSided;
    pixelsHalfAngleDoubleSided.resize(IMAGESIZE()*IMAGESIZE() * 4);
    std::fill(pixelsHalfAngleDoubleSided.begin(), pixelsHalfAngleDoubleSided.end(), 255);
    DrawCircle(pixelsHalfAngleDoubleSided, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 2, 240, 240, 240);
    DrawCircle(pixelsHalfAngleDoubleSided, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, IMAGESIZE() / 2 - 4, 255, 255, 255);

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

            DrawLine(pixelsAngle, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, lineEndX, lineEndY, color, color, color);

            // also draw the double sided version
            DrawLine(pixelsAngleDoubleSided, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() - lineEndX, IMAGESIZE() - lineEndY, lineEndX, lineEndY, color, color, color);
        }

        // use golden ratio to make angles between 0 and pi
        {
            float angleRad = std::fmodf(float(i) * c_goldenRatioConjugate, 1.0f) * c_pi;
            float angleDeg = angleRad * 180.0f / c_pi;

            int lineEndX = int(cos(angleRad) * (float(IMAGESIZE() / 2) - 4) + float(IMAGESIZE() / 2));
            int lineEndY = int(sin(angleRad) * (float(IMAGESIZE() / 2) - 4) + float(IMAGESIZE() / 2));

            lineEndY = IMAGESIZE() - lineEndY;

            uint8 color = uint8(255.0f * float(i) / float(c_numAngles));

            DrawLine(pixelsHalfAngle, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() / 2, IMAGESIZE() / 2, lineEndX, lineEndY, color, color, color);

            // also draw the double sided version
            DrawLine(pixelsHalfAngleDoubleSided, IMAGESIZE(), IMAGESIZE(), IMAGESIZE() - lineEndX, IMAGESIZE() - lineEndY, lineEndX, lineEndY, color, color, color);
        }
    }

    SaveImage("out/AngleTest.png", IMAGESIZE(), IMAGESIZE(), pixelsAngle.data());
    SaveImage("out/AngleTestDS.png", IMAGESIZE(), IMAGESIZE(), pixelsAngleDoubleSided.data());
    SaveImage("out/AngleTestHalf.png", IMAGESIZE(), IMAGESIZE(), pixelsHalfAngle.data());
    SaveImage("out/AngleTestHalfDS.png", IMAGESIZE(), IMAGESIZE(), pixelsHalfAngleDoubleSided.data());
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
    LHS lhs_blueNoise(points, GeneratePointHashDatas_BlueNoise);
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
        ReportQuery(lhs_blueNoise, queryPoint, points, "blueNoise");
        ReportQuery(lhs_uniform, queryPoint, points, "uniform");
        ReportQuery(lhs_goldenRatio, queryPoint, points, "goldenRatio");
        ReportQuery(lhs_goldenRatioGeneralized, queryPoint, points, "goldenRatioGeneralized");
    }

    // TODO: ground truth, how? maybe visually show points missed or something? could calculate std deviations as concentric rings.

    return 0;
}

/*
 TODO:
 * tune blue noise now that you can see it plotted.
 * 3 blue noises? 1) blue noise independant on each axis.  2) 2d blue noise. 3) projective blue noise on x and y axis. 4) ??? random axis projective blue noise?
* uniform needs help. also i think it is incorrectly calculating how many rows it needs. set num hashes to 20 and note how it doesn't look correct.
* make a light image class that stores pixels, width, height
 * it would be neat to make a header only library for this image class and primitive drawing, for use in other blog posts too.
* test the offset specifically, with 1 bucket and an offset. make sure it's happy :)
? are the buckets too small? i think they might be... instead of having the domain thing for points, could have a bucket size scale.
* maybe an option to color cells based on how many hash collisions they have with the query point
* tyler's article says there is an optimal for white noise. check out what that is, and try it? maybe compare things with that optimal value?
* make a function to save png's, so you can put your "write a transparent pixel" code in one place :p
* there is likely a point where white noise beats out blue / LDS. find it?

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
 * the paper from tyler sort of mentions a bloom filter related thing: "if we’re willing to probabilistically verify hash keys at the cost of some false collisions, we can still achieve O(d log d) lookups" 
 * basically, if you get hashed items, you need to resolve them to points which takes longer. if you work only w/ hashes you can get false positives but it's faster. A bit bloom filter ish.

* there was also a technique talked about


Links:
- inspiration: http://unboxresearch.com/articles/lsh_post1.html   aka http://unboxresearch.com/articles/lsh_post1.html
* another way to do it: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.215.7690&rep=rep1&type=pdf
 * make regular cells of any shape. store a list at each cell corner which has each point that is in a cell that uses that corner.
 * talks about simplices fighting the curse of dimensionality.
 * regular simplices can't tile higher dimensions though so it has an alternate formulation.
- maybe some blue noise articles?
* frobenius inner product for checking similarity of rotation matrices: https://en.wikipedia.org/wiki/Frobenius_inner_product
* random rotation matrices: https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices

Notes:
* have a white noise [0,1] random value to start golden ratio by. You get the benefits of random values i believe, like you do with white noise.
 * do you get all of them?
* Does it make sense to do golden ratio between 0 and pi instead of 0 and 2 pi? yes. check the image tests. The angles are "double sided".
* mitchell's best candidate algorithm for generating blue noise has a tuneable parameter. Too low and you get white noise. Too high and you get regular sampling. not super great for getting good results out of the box.
* generalized golden ratio doesn't look that LDS at lower hash counts, but it does at higher hash counts.  Maybe not useful since hash count should (?) always be same as # of faces on a simplex for that dimension?

----- LANDFILL -----

* do we need the full matrix? i don't think so, if we are only takiung the resulting x component!
 * maybe we don't care. implementation detail (optimization)

*/