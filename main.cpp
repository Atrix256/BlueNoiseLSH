// Configurable settings
#define DIMENSION() 2      // dimensionality of the data points
#define NUMPOINTS() 100    // the number of randomly generated (white noise) data points
#define HASHCOUNT() 4      // how many hashes are done for each point


#include <random>
#include <array>
#include <vector>
#include <unordered_set>
#include <stdint.h>

typedef std::array<float, DIMENSION()> TPoint;
typedef uint32_t uint32;

enum class PointID : uint32 {};

std::mt19937& GetRNG()
{
    static std::random_device rd;
    static std::mt19937 rng(rd());
    return rng;
}

class LHS
{
public:
    LHS (const std::vector<TPoint>& points)
    {
        m_points = points;
        for (const TPoint& point : m_points)
            AddPoint(point);
    }

private:

    void AddPoint (const TPoint& point)
    {

    }

private:

    std::vector<TPoint> m_points; // a point ID (uint32) is the point's index in this list

    std::array<std::unordered_set<PointID>, HASHCOUNT()> m_hashBuckets;
};

int main(int argc, char** argv)
{
    std::mt19937& rng = GetRNG();
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

    // make the random point set
    std::vector<TPoint> points;
    points.resize(NUMPOINTS());
    for (TPoint& p : points)
    {
        for (float& f : p)
            f = dist(rng);
    }

    // add the points to the locality sensitive hashing
    LHS lhs(points);



    return 0;
}

/*

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

*/