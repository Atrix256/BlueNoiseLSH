// Configurable settings
#define DIMENSION() 2      // dimensionality of the data points
#define NUMPOINTS() 100    // the number of randomly generated (white noise) data points
#define HASHCOUNT() 4      // how many hashes are done for each point
#define POINTDOMAIN() 10   // the coordinates go from - this value to + this value

#include <random>
#include <array>
#include <vector>
#include <unordered_set>
#include <stdint.h>

// -------------------------------------------------------------------------------

typedef std::array<float, DIMENSION()> TPoint;
typedef uint32_t uint32;

enum class PointID : uint32 {};

static const float c_pi = 3.14159265359f;

// -------------------------------------------------------------------------------

struct PointHashData
{
    std::array<float, DIMENSION()*DIMENSION()> rotation;
    float offsetX;
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

private:

    void AddPoint (const TPoint& point, PointID pointID)
    {
        for (int hashIndex = 0; hashIndex < HASHCOUNT(); ++hashIndex)
        {
            const PointHashData& pointHashData = m_pointHashDatas[hashIndex];

            TPoint rotatedPoint;

            for (int i = 0; i < DIMENSION(); ++i)
            {
                int rowOffset = i * DIMENSION();

                rotatedPoint[i] = 0.0f;
                for (int j = 0; j < DIMENSION(); ++j)
                    rotatedPoint[i] += point[j] * pointHashData.rotation[rowOffset + j];
            }

            int hashValue = (int)std::floorf(rotatedPoint[0] + pointHashData.offsetX);

            // TODO: put the point id into the bucket defined by hashValue - may be negative!
            int ijkl = 0;
        }
    }

private:

    std::vector<TPoint> m_points; // a PointID (uint32) is the point's index in this list

    std::array<PointHashData, HASHCOUNT()> m_pointHashDatas;

    //std::array<std::unordered_set<PointID>, HASHCOUNT()> m_hashBuckets;
};

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
        float sinTheta = std::cosf(angle);

        p.rotation = { cosTheta, -sinTheta,
                       sinTheta,  cosTheta };

        p.offsetX = dist_offset(RNG());
    }
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
    LHS lhs(points, GeneratePointHashDatas_WhiteNoise);



    return 0;
}

/*

* do we need the full matrix? i don't think so, if we are only takiung the resulting x component!

Tests:
* white noise
* uniform
* blue noise
* LDS (golden ratio?)
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