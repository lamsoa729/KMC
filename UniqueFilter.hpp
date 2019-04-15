/**
 * @file UniqueFilter.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-02-06
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef UNIQUEFILTER_HPP_
#define UNIQUEFILTER_HPP_

#include <unordered_set>
#include <vector>

/**
 * @brief uniqueness filter for FDPS interaction list.
 * objs with the same gid AND the same position are considered the same object.
 *
 * @tparam T
 */
template <class T>
class UniqueFilter {
  private:
    const T *const objPtr; ///< pointer to obj data
    const int N;        ///< number
    std::vector<int> uniqueFlag;

    /**
     * @brief Entry in unordered map
     *
     */
    struct Entry {
        int gid;
        double x, y, z;

        explicit Entry(const T &obj) {
            gid = obj.getGid();
            const auto &vec = obj.getPos();
            x = vec.x;
            y = vec.y;
            z = vec.z;
        }
    };

    /**
     * @brief Combined hash functor
     *
     */
    class HashFN {
      public:
        size_t operator()(const Entry &obj) const {
            size_t hash = 0;
            std::hash<int> intHasher;
            std::hash<double> doubleHasher;
            hash = intHasher(obj.gid) * 0.25 + doubleHasher(obj.x) * 0.25 +
                   doubleHasher(obj.y) * 0.25 + doubleHasher(obj.z) * 0.25;
            return hash;
        }
    };

    /**
     * @brief Equality functor
     *
     */
    class EqualFN {
      public:
        bool operator()(const Entry &obj1, const Entry &obj2) const {
            auto rx = obj1.x - obj2.x;
            auto ry = obj1.y - obj2.y;
            auto rz = obj1.z - obj2.z;
            if (obj1.gid == obj2.gid && rx * rx + ry * ry + rz * rz < 1e-8) {
                return true;
            } else {
                return false;
            }
        }
    };

    /**
     * @brief update the unique flag for obj data
     *
     * @param obj
     * @param Nobj
     */
    void updateUniqueFlag() {
        if (N <= 0) {
            uniqueFlag.clear();
            return;
        }

        std::unordered_set<Entry, HashFN, EqualFN> gidSet;
        uniqueFlag.resize(N, 0);

        for (int i = 0; i < N; i++) {
            auto &objI = objPtr[i];
            Entry entry(objI);
            if (gidSet.find(entry) != gidSet.end()) {
                uniqueFlag[i] = -1; // not unique, repeated
            } else {
                gidSet.insert(entry);
                uniqueFlag[i] = 1; // unique or first of non unique
            }
        }
    }

  public:
    /**
     * @brief Construct a new UniqueFilter object
     *
     * @param objPtr_
     * @param N_
     */
    UniqueFilter(const T *const objPtr_, int N_) : objPtr(objPtr_), N(N_) {
        uniqueFlag.resize(N);
        updateUniqueFlag();
    }

    ~UniqueFilter() = default;

    /**
     * @brief Get the Filter object
     *
     * @return std::vector<int>&
     */
    std::vector<int> &getFilter() { return uniqueFlag; }
};

#endif
