#include "FixedKmerGenerator.h"
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_set>
#include <numeric>
#include <tuple>

FixedKmerGenerator::FixedKmerGenerator(size_t kmerSize, size_t alphabetSize, short threshold, unsigned int maxKmers)
    : kmerSize(kmerSize), threshold(threshold), maxKmers(maxKmers), indexer((int) alphabetSize, (int)kmerSize) {}

FixedKmerGenerator::~FixedKmerGenerator() {
    delete[] stepMultiplicator;
    delete[] divideStep;
    delete[] matrixLookup;
    delete[] outputIndexArray;
}

void FixedKmerGenerator::setThreshold(short threshold) {
    this->threshold = threshold;
}

void FixedKmerGenerator::setDivideStrategy(ScoreMatrix** one) {
    divideStepCount = kmerSize;
    matrixLookup = new ScoreMatrix*[divideStepCount];
    divideStep   = new unsigned int[divideStepCount];
    for (size_t i = 0; i < kmerSize; i++) {
        divideStep[i] = 1;
        matrixLookup[i] = one[i];
    }
    initDataStructure();
}

void FixedKmerGenerator::setDivideStrategy(ScoreMatrix* three, ScoreMatrix* two) {
    const size_t threeDivideCount = kmerSize / 3;
    switch (kmerSize % 3) {
        case 0:
            divideStepCount = threeDivideCount;
            matrixLookup = new ScoreMatrix*[divideStepCount];
            divideStep = new unsigned int[divideStepCount];
            for (size_t i = 0; i < threeDivideCount; i++) {
                divideStep[i] = 3;
                matrixLookup[i] = three;
            }
            break;
        case 1:
            divideStepCount = threeDivideCount+1;
            matrixLookup = new ScoreMatrix*[divideStepCount];
            divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount-1; i++){
                divideStep[i] = 3;
                matrixLookup[i] = three;
            }
            divideStep[threeDivideCount-1] = 2;
            matrixLookup[threeDivideCount-1] = two;

            divideStep[threeDivideCount] = 2;
            matrixLookup[threeDivideCount] = two;

            break;
        case 2:
            divideStepCount = threeDivideCount+1;
            matrixLookup = new ScoreMatrix*[divideStepCount];
            divideStep = new unsigned int[divideStepCount];
            for(size_t i = 0; i < threeDivideCount; i++){
                divideStep[i] = 3;
                matrixLookup[i] = three;
            }
            divideStep[threeDivideCount] = 2;
            matrixLookup[threeDivideCount] = two;

            break;
    }

    initDataStructure();
    std::reverse(matrixLookup, &matrixLookup[divideStepCount]);
    std::reverse(divideStep, &divideStep[divideStepCount]);
}


void FixedKmerGenerator::initDataStructure() {
    stepMultiplicator = new size_t[divideStepCount];
    outputIndexArray = (size_t *) mem_align(ALIGN_INT, maxKmers * sizeof(size_t));
}

// std::pair<std::vector<std::vector<short>>, std::vector<std::vector<unsigned int>>> topScoringKmers(std::vector<std::pair<short*, int>> &scores, std::vector<std::pair<unsigned int*, int>> &indicesPairs, int N) {
//     // std::priority_queue<std::pair<int, std::vector<int>>, std::vector<std::pair<int, std::vector<int>>>, decltype(&comp)> queue(&comp);
//     std::vector<std::pair<int, std::vector<int>>> queue;
//     std::make_heap(queue.begin(), queue.end(), comp);
//     std::set<std::vector<int>> visited;

//     std::vector<int> indices(scores.size(), 0);

//     int total = 0;
//     for (size_t i = 0; i < scores.size(); ++i) {
//         total += scores[i].first[0];
//     }

//     // queue.emplace(total, indices);
//     queue.push_back(std::make_pair(total, indices));
//     std::push_heap(queue.begin(), queue.end(), comp);
//     visited.insert(indices);

//     std::vector<std::vector<short>> result;
//     std::vector<std::vector<unsigned int>> indexResult;

//     for (int _ = 0; _ < N; ++_) {
//         // const std::pair<int, std::vector<int>>& top = queue.top();
//         // total = top.first;
//         // indices = top.second;
//         // queue.pop();
//         std::pop_heap(queue.begin(), queue.end(), comp);
//         std::pair<int, std::vector<int>>& back = queue.back();
//         total = back.first;
//         indices = back.second;
//         queue.pop_back();

//         std::vector<short> res;
//         std::vector<unsigned int> indexRes;
//         for (size_t i = 0; i < scores.size(); ++i) {
//             res.push_back(scores[i].first[indices[i]]);
//             indexRes.push_back(indicesPairs[i].first[indices[i]]);
//         }
//         result.push_back(res);
//         indexResult.push_back(indexRes);

//         for (size_t i = 0; i < scores.size(); ++i) {
//             if (indices[i] + 1 < scores[i].second) {
//                 std::vector<int> newIndices = indices;
//                 newIndices[i] += 1;
//                 if (visited.find(newIndices) == visited.end()) {
//                     visited.insert(newIndices);
//                     total = 0;
//                     for (size_t j = 0; j < scores.size(); ++j) {
//                         total += scores[j].first[newIndices[j]];
//                     }
//                     // queue.emplace(total, newIndices);
//                     queue.push_back(std::make_pair(total, newIndices));
//                     std::push_heap(queue.begin(), queue.end(), comp);
//                 }
//             }
//         }
//     }

//     return {result, indexResult};
// }

template <typename T>
bool comp(const T& a, const T& b) {
    return std::get<0>(a) < std::get<0>(b);
}

struct pair_hash_fn {
    template <class T>
    inline void hash_combine(std::size_t& seed, T const& v) const {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    template <class Pair>
    inline std::size_t operator()(Pair const& t) const
    {
        std::size_t seed = 0;
        hash_combine(seed, t.first);
        hash_combine(seed, t.second);
        return seed;
    }
};

struct tuple_hash_fn {
    template <class T>
    inline void hash_combine(std::size_t& seed, T const& v) const
    {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    template <class Tuple>
    inline std::size_t operator()(Tuple const& t) const {
        std::size_t seed = 0;
        hash_combine(seed, std::get<0>(t));
        hash_combine(seed, std::get<1>(t));
        hash_combine(seed, std::get<2>(t));
        return seed;
    }
};


std::pair<std::vector<unsigned int>, std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>> topScoringKmers2(std::pair<short*, int>* scores,std::pair<unsigned int*, int>* indices, int N) {
    std::vector<std::tuple<int, int, int>> queue;
    std::unordered_set<std::pair<int, int>, pair_hash_fn> visited;

    queue.emplace_back(scores[0].first[0] + scores[1].first[0], 0, 0);
    std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int>>);
    visited.emplace(0, 0);

    std::vector<unsigned int> result;
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> indexResult;

    for (int _ = 0; _ < N; ++_) {
        std::pop_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int>>);
        const std::tuple<int, int, int>& back = queue.back();
        int total = std::get<0>(back);
        int i = std::get<1>(back);
        int j = std::get<2>(back);
        queue.pop_back();

        result.emplace_back(total);
        indexResult.emplace_back(indices[0].first[i], indices[1].first[j], 0);

        if (i + 1 < scores[0].second && visited.find({i + 1, j}) == visited.end()) {
            queue.emplace_back(scores[0].first[i + 1] + scores[1].first[j], i + 1, j);
            std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int>>);
            visited.emplace(i + 1, j);
        }

        if (j + 1 < scores[1].second && visited.find({i, j + 1}) == visited.end()) {
            queue.emplace_back(scores[0].first[i] + scores[1].first[j + 1], i, j + 1);
            std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int>>);
            visited.emplace(i, j + 1);
        }
    }

    return {result, indexResult};
}

std::pair<std::vector<unsigned int>, std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>> topScoringKmers3(std::pair<short*, int>* scores, std::pair<unsigned int*, int>* indices, int N) {
    std::vector<std::tuple<int, int, int, int>> queue;
    std::unordered_set<std::tuple<int, int, int>, tuple_hash_fn> visited;

    queue.emplace_back(scores[0].first[0] + scores[1].first[0] + scores[2].first[0], 0, 0, 0);
    std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int, int>>);
    visited.emplace(0, 0, 0);

    std::vector<unsigned int> result;
    result.reserve(N);
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> indexResult;
    indexResult.reserve(N);

    for (int _ = 0; _ < N; ++_) {
        std::pop_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int, int>>);
        const std::tuple<int, int, int, int>& back = queue.back();
        int total = std::get<0>(back);
        int i = std::get<1>(back);
        int j = std::get<2>(back);
        int k = std::get<3>(back);
        queue.pop_back();

        result.emplace_back(total);
        indexResult.emplace_back(indices[0].first[i], indices[1].first[j], indices[2].first[k]);

        if (i + 1 < scores[0].second && visited.find({i + 1, j, k}) == visited.end()) {
            queue.emplace_back(scores[0].first[i + 1] + scores[1].first[j] + scores[2].first[k], i + 1, j, k);
            std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int, int>>);
            visited.emplace(i + 1, j, k);
        }

        if (j + 1 < scores[1].second && visited.find({i, j + 1, k}) == visited.end()) {
            queue.emplace_back(scores[0].first[i] + scores[1].first[j + 1] + scores[2].first[k], i, j + 1, k);
            std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int, int>>);
            visited.emplace(i, j + 1, k);
        }

        if (k + 1 < scores[2].second && visited.find({i, j, k + 1}) == visited.end()) {
            queue.emplace_back(scores[0].first[i] + scores[1].first[j] + scores[2].first[k + 1], i, j, k + 1);
            std::push_heap(queue.begin(), queue.end(), comp<std::tuple<int, int, int, int>>);
            visited.emplace(i, j, k + 1);
        }
    }

    return {result, indexResult};
}

std::pair<size_t *, size_t> FixedKmerGenerator::generateKmerList(const unsigned char * int_seq, bool /* addIdentity */) {
    int dividerBefore = 0;
    std::pair<short *, int> scoreArrays[3];
    std::pair<unsigned int*, int> indexArrays[3];
    for (size_t i = 0; i < divideStepCount; i++) {
        const int divider = divideStep[i];
        const unsigned int index = indexer.int2index(int_seq, dividerBefore, dividerBefore+divider);
        stepMultiplicator[i] = indexer.powers[dividerBefore];
        dividerBefore += divider;

        const ScoreMatrix* nextScoreMatrix = matrixLookup[i];
        short* nextScoreArray = &nextScoreMatrix->score[index*nextScoreMatrix->rowSize];
        scoreArrays[i] = std::make_pair(nextScoreArray, (int)nextScoreMatrix->rowSize);

        unsigned int* nextIndexArray = &nextScoreMatrix->index[index*nextScoreMatrix->rowSize];
        indexArrays[i] = std::make_pair(nextIndexArray, (int)nextScoreMatrix->rowSize);
    }
    
    std::pair<std::vector<unsigned int>, std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>> res;
    if (divideStepCount == 2) {
        res = topScoringKmers2(scoreArrays, indexArrays, maxKmers);
    } else if (divideStepCount == 3)
        res = topScoringKmers3(scoreArrays, indexArrays, maxKmers);
    else {
        std::cerr << "Unsupported divideStepCount: " << divideStepCount << std::endl;
        exit(1);
    }

    size_t generated = 0;
    for (size_t i = 0; i < res.first.size(); i++) {
        unsigned int score = res.first[i];
        size_t index = 0;
        index += static_cast<size_t>(std::get<0>(res.second[i])) * stepMultiplicator[0];
        index += static_cast<size_t>(std::get<1>(res.second[i])) * stepMultiplicator[1];
        if (divideStepCount == 3) {
            index += static_cast<size_t>(std::get<2>(res.second[i])) * stepMultiplicator[2];
        }
        if (score >= (unsigned int)threshold) {
            outputIndexArray[generated++] = index;
        }
    }

    return std::make_pair(outputIndexArray, generated);
}
