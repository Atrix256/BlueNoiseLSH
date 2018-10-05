// Configurable settings
#define DIMENSION() 2      // dimensionality of the data points
#define NUMPOINTS() 100    // the number of randomly generated (white noise) data points
#define HASHCOUNT() 1      // how many hashes are done for each point
#define POINTDOMAIN() 10   // the coordinates go from - this value to + this value

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

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = 0.0f;// dist_offset(RNG());
    }
}

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

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = 0.0f;// percentY;
    }
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

            pixel[0] = Lerp(pixel[0], R, alpha);
            pixel[1] = Lerp(pixel[1], G, alpha);
            pixel[2] = Lerp(pixel[2], B, alpha);

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
        for (int ix = 0; ix < imageWidth; ++ix)
        {
            // transform pixel into the hash bucket space
            TPoint rawPoint = {float(ix), float(iy)};
            TPoint rotatedPoint = Multiply(rawPoint, pointHashData.rotation);
            rotatedPoint[0] += pointHashData.offsetX;

            // figure out the distance from this rotated point to the edges of the hash bucket
            float cellPosX;
            if (rotatedPoint[0] < 0.0f)
                cellPosX = cellSize - fmodf(fabsf(rotatedPoint[0]) + 0.5f, cellSize);
            else
                cellPosX = fmodf(rotatedPoint[0] + 0.5f, cellSize);
            float distance = cellSize / 2.0f - fabsf(cellPosX - cellSize / 2.0f);

            // draw the hash bucket lines
            float alpha = SmoothStep(distance, 2.0f, 0.0f);
            pixel[0] = Lerp(pixel[0], R, alpha);
            pixel[1] = Lerp(pixel[1], G, alpha);
            pixel[2] = Lerp(pixel[2], B, alpha);
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

    int width = 500;
    int height = 500;
    int circleRadius = int(float(width) / 100.0f);
    int colorCircleRadius = std::min(int(float(width) / 100.0f * 0.9f), circleRadius - 1);

    // create the image and fill it with white
    std::vector<uint8> pixels;
    pixels.resize(width*height * 4);
    std::fill(pixels.begin(), pixels.end(), 255);

    // draw the hash buckets
    for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
    {
        const PointHashData& pointHashData = lhs.GetPointHashData(hashIndex);
        DrawHashBuckets(pixels, width, height, float(width) * 0.5f / float(POINTDOMAIN()), pointHashData, 128, 128, 128);
    }

    // draw all the points as black dots
    for (TPoint p : points)
    {
        for (float& f : p)
        {
            f /= (2.0f * float(POINTDOMAIN()));
            f += 0.5f;
        }

        int x = int(p[0] * float(width - 1) + 0.5f);
        int y = int(p[1] * float(height - 1) + 0.5f);

        DrawCircle(pixels, width, height, x, y, circleRadius, 0, 0, 0);
    }

    // draw the results based on their match count
    for (auto it : results)
    {
        TPoint p = points[(uint32)it.first];

        for (float& f : p)
        {
            f /= (2.0f * float(POINTDOMAIN()));
            f += 0.5f;
        }

        int x = int(p[0] * float(width - 1) + 0.5f);
        int y = int(p[1] * float(height - 1) + 0.5f);

        float percentMatched = float(it.second) / float(HASHCOUNT());
        uint8 color = uint8(percentMatched * 255.0f);

        DrawCircle(pixels, width, height, x, y, colorCircleRadius, 255, color, 0);
    }

    // draw the query point
    {
        TPoint p = queryPoint;
        for (float& f : p)
        {
            f /= (2.0f * float(POINTDOMAIN()));
            f += 0.5f;
        }
        int x = int(p[0] * float(width - 1) + 0.5f);
        int y = int(p[1] * float(height - 1) + 0.5f);

        DrawCircle(pixels, width, height, x, y, circleRadius, 0, 0, 255);
    }

    stbi_write_png(fileName, width, height, 4, pixels.data(), 0);
}

// -------------------------------------------------------------------------------

void ReportQuery (const LHS& lhs, const TPoint& queryPoint, const std::vector<TPoint>& points, const char* name)
{
    printf("\n=====%s=====\n", name);

    std::unordered_map<PointID, int> results;
    lhs.Query(queryPoint, results);

    for (int hashMatchCount = HASHCOUNT(); hashMatchCount >= 1; --hashMatchCount)
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

int main(int argc, char** argv)
{
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
    LHS lhs_white(points, GeneratePointHashDatas_WhiteNoise);
    LHS lhs_uniform(points, GeneratePointHashDatas_Uniform);

    // do some queries
    for (int i = 0; i < 1; ++i)
    {
        TPoint queryPoint;
        for (float& f : queryPoint)
            f = dist(RNG());

        printfPoint("Query Point", queryPoint);

        ReportQuery(lhs_white, queryPoint, points, "white");
        ReportQuery(lhs_uniform, queryPoint, points, "uniform");
    }

    // TODO: ground truth, how? maybe visually show points missed or something? could calculate std deviations as concentric rings.

    return 0;
}

// TODO: it seems like the grid lines may be wrong sometimes even with no x offset. check out the white noise screenshot checked in. It probably has to do with where it goes negative. the distance calculation is probably wrong...

/*

TODO:
* make the offset work with grid drawing. I think it has something to do with how we floor it when we use it for points so it's a randomized rounding? dunno. also seems to maybe happen when the point is a negative number?
* need a minimum hash count i think. make a #define
* maybe an option to color cells based on how many hash collisions they have with the query point
? is the grid line at 0 thinner? i think it is... why?

* do we need the full matrix? i don't think so, if we are only takiung the resulting x component!
 * maybe we don't care

Tests:
* white noise
* uniform
* blue noise
* LDS (golden ratio?) - yes. generalized golden ratio. http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
* compare vs linear search to get ground truth

* Some way to gather stats, but also want to make diagrams i think. Maybe some animated like the other blog post does?

maybe MATCHHASHCOUNT should be a parameter to a point query function

----- LATER -----


* maybe a second blog post on bloom filter?
 * blue noise is about quality of results.
 * bloom filter is about performance, but may also affect quality
 * the paper from tyler sort of mentions a bloom filter related thing: "if we’re willing to probabilistically verify hash keys at the cost of some false collisions, we can still achieve O(d log d) lookups" 
 * basically, if you get hashed items, you need to resolve them to points which takes longer. if you work only w/ hashes you can get false positives but it's faster. A bit bloom filter ish.

* there was also a technique talked about


Links:
- inspiration: http://unboxresearch.com/articles/lsh_post1.html
* another way to do it: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.215.7690&rep=rep1&type=pdf
 * make regular cells of any shape. store a list at each cell corner which has each point that is in a cell that uses that corner.
 * talks about simplices fighting the curse of dimensionality.
 * regular simplices can't tile higher dimensions though so it has an alternate formulation.
- maybe some blue noise articles?
* random rotation matrices: https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices

*/