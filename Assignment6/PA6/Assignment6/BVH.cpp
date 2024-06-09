#include <algorithm>
#include <array>
#include <cassert>
#include <cfloat>
#include <vector>
#include "BVH.hpp"
#include "Bounds3.hpp"
#include "Object.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else if (splitMethod == SplitMethod::NAIVE) {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    // TODO: implement SAH accelerating structure
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        // TODO: Patitioning by buckets
        float SN = centroidBounds.SurfaceArea();
        double leftside, rightside;
        switch (dim) {
            case 0:
                leftside = centroidBounds.pMin.x;
                rightside = centroidBounds.pMax.x;
                break;
            case 1:
                leftside = centroidBounds.pMin.y;
                rightside = centroidBounds.pMax.y;
                break;
            case 2:
                leftside = centroidBounds.pMin.z;
                rightside = centroidBounds.pMax.z;
                break;
        }
        int B = 16;
        auto step = (rightside - leftside) / (float)B;
        int bucketCount[B];
        int minCostIndex = 0;
        float minCost = FLT_MAX;
        Bounds3 leftBounds, rightBounds;
        for (int i = 1; i < B; i++)
        {
            auto beginning = objects.begin();
            auto ending = objects.end();
            auto middling = objects.begin();
            int index = 0;
            if (i > 1)
            {
                middling += bucketCount[i];
                index = bucketCount[i];
            }

            while (middling != ending && dim == 0
                       ? ((*(middling + 1))->getBounds().Centroid()).x
                       : (dim == 1
                              ? ((*(middling + 1))->getBounds().Centroid()).y
                              : ((*(middling + 1))->getBounds().Centroid()).z)
                             < leftside + i * step)
            {
                middling++;
                index++;
            }
            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);
            assert(objects.size() == (leftshapes.size() + rightshapes.size()));
            for (int k = 0; k < leftshapes.size(); k++)
                leftBounds
                    = Union(leftBounds, leftshapes[k]->getBounds().Centroid());
            for (int k = 0; k < rightshapes.size(); k++)
                rightBounds = Union(rightBounds,
                                    rightshapes[k]->getBounds().Centroid());
            float SA = leftBounds.SurfaceArea();
            float SB = rightBounds.SurfaceArea();
            float cost
                = 0.125
                  + (leftshapes.size() * SA + rightshapes.size() * SB) / SN;
            if (cost < minCost)
            {
                minCost = cost;
                minCostIndex = index;
            }
        }

        auto beginning = objects.begin();
        auto middling = beginning + minCostIndex;
        auto ending = objects.end();
        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));
        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);
        node->bounds = Union(node->left->bounds, node->right->bounds);
    }
    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    std::array<int, 3> dirIsNeg;
    dirIsNeg[0] = ray.direction.x > 0 ? 0 : 1;
    dirIsNeg[1] = ray.direction.y > 0 ? 0 : 1;
    dirIsNeg[2] = ray.direction.z > 0 ? 0 : 1;

    if (!node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg))
        return {};

    if (node->left==nullptr && node->right==nullptr) {
        return node->object->getIntersection(ray);
    }
    
    Intersection leaf1 = BVHAccel::getIntersection(node->left, ray);
    Intersection leaf2 = BVHAccel::getIntersection(node->right, ray);
    return leaf1.distance < leaf2.distance ? leaf1 : leaf2;
}