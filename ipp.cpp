/**
 * C++ implementation of the ipp module.
 */
#include "ipp.h"

// Always execute assert()s!
#undef NDEBUG

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <tuple>
#include <unordered_set>

namespace {

template<typename T>
T
readInt(std::ifstream& file) {
    T val;
    file.read(reinterpret_cast<char*>(&val), sizeof(val));
    if (!file.good()) {
        throw std::runtime_error("Unexpected EOF");
    }
    return val;
}

std::string
readString(std::ifstream& file) {
    std::string s;
    std::getline(file, s, '\0');
    if (!file.good()) {
        throw std::runtime_error("Unexpected EOF");
    }
    return s;
}

} // namespace

void
Ipp::loadPwalns(std::string const& fileName) {
    // Reads the chromosomes and pwalns from the given file.
    // The data in the file is expected to be in the following format:
    //
    // version                   [uint8]
    // endianness_magic          [uint16]
    // num_sp1                   [uint8]
    // {
    //   sp1_name                [null-terminated string]
    //   sp1_genome_size         [uint64]
    //   num_sp2                 [uint8]
    //   {
    //     sp2_name              [null-terminated string]
    //     num_ref_chrom_entries [uint32]
    //     {
    //       ref_chrom           [uint32]
    //       num_pwaln_entries   [uint32]
    //       {
    //         ref_start         [uint32]
    //         qry_start         [uint32]
    //         qry_chrom         [uint32]
    //         length_and_strand [uint16]
    //       } num_pwaln_entries times
    //     } num_ref_chrom_entries times
    //   } num_sp2 times
    // } num_sp1 times
    // num_chomosomes            [uint32]
    // {
    //   chrom_name              [null-terminated string]
    // } num_chromosomes times

    uint8_t expectedFormatVersion(4);

    chroms_.clear();
    genomeSizes_.clear();
    pwalns_.clear();

    std::ifstream file(fileName, std::ios::in|std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("could not open the file");
    }

    // Read and check the format version.
    auto const formatVersion(readInt<uint8_t>(file));
    if (formatVersion != expectedFormatVersion) {
        throw std::runtime_error(format("invalid version: %u (expected: %u)",
                                        formatVersion, expectedFormatVersion));
    }

    // Read and check the endianness magic number.
    auto const endiannessMagic(readInt<uint16_t>(file));
    if (endiannessMagic != 0xAFFE) {
        throw std::runtime_error(
            "the endianness of the system that produced the pwalns file "
            "differs from the enndianess of this system");
    }

    // Read the pwalns.
    maxAnchorLength_ = 0;
    auto const numSp1(readInt<uint8_t>(file));
    for (unsigned i(0); i < numSp1; ++i) {
        // Create new entry in the map if it is not yet there.
        std::string const sp1(readString(file));

        genomeSizes_[sp1] = readInt<uint64_t>(file);
        auto& pwalnsSp1(pwalns_[sp1]);

        auto const numSp2(readInt<uint8_t>(file));
        for (unsigned j(0); j < numSp2; ++j) {
            std::string const sp2(readString(file));
            auto& pwaln(pwalnsSp1[sp2]);

            auto const numRefChromEntries(readInt<uint32_t>(file));
            for (unsigned k(0); k < numRefChromEntries; ++k) {
                auto const refChrom(readInt<uint32_t>(file));
                auto const numPwalnEntries(readInt<uint32_t>(file));

                // Bulk-read the pwaln entries into a temporary buffer.
                class BufAnchor {
                public:
                   explicit BufAnchor(std::size_t size)
                       : size(size)
                       , buf(std::malloc(size))
                   {}
                   ~BufAnchor() {
                       std::free(buf);
                   }

                   std::size_t const size;
                   void* const buf;
                };

                #pragma pack(1)
                class PackedPwalnEntry : public PwalnEntry {
                    // Packed PwalnEntry that does not have the extra padding
                    // bytes at the end which are used for alignment.
                };
                static_assert(sizeof(PackedPwalnEntry) == 14);
                #pragma pack()

                BufAnchor bufAnchor(numPwalnEntries*sizeof(PackedPwalnEntry));
                auto const buf(
                    reinterpret_cast<PackedPwalnEntry*>(bufAnchor.buf));

                file.read(reinterpret_cast<char*>(buf), bufAnchor.size);
                if (!file.good()) {
                    throw std::runtime_error("Unexpected EOF");
                }

                // Copy the packed representation to the properly aligned
                // vector.
                static_assert(sizeof(PwalnEntry) == 16);
                std::vector<PwalnEntry>& pwalnEntries(pwaln[refChrom]);
                pwalnEntries.reserve(numPwalnEntries);
                for (unsigned l(0); l < numPwalnEntries; ++l) {
                    maxAnchorLength_ = std::max(maxAnchorLength_,
                                                buf[l].length());
                    pwalnEntries.emplace_back(buf[l]);
                }
            }
        }
    }

	// Read the chromosomes.
    auto const numChromosomes(readInt<uint32_t>(file));
    chroms_.reserve(numChromosomes);
    for (unsigned i(0); i < numChromosomes; ++i) {
        chroms_.push_back(readString(file));
    }

    if (file.peek() != std::ifstream::traits_type::eof()) {
        throw std::runtime_error("Remaining data at EOF");
        // There is more data to read when we don't expect it.
    }
}

void
Ipp::setHalfLifeDistance(unsigned halfLifeDistance) {
    // Sets the half-life distance;
    halfLifeDistance_ = halfLifeDistance;
}

std::optional<Ipp::ChromId>
Ipp::chromIdFromName(std::string const& chromName) const {
    // Looks up the given chromosome name in chroms_ and returns its id.
    auto const it(std::find(chroms_.begin(), chroms_.end(), chromName));
    if (it == chroms_.end()) {
        return {};
    }
    return {std::distance(chroms_.begin(), it)};
}

std::string const&
Ipp::chromName(ChromId chromId) const {
    // Returns the name of the chromosome with the given id.
    return chroms_.at(chromId);
}

double
Ipp::projectionScore(uint32_t loc,
                     uint32_t upBound,
                     uint32_t downBound,
                     uint64_t genomeSize,
                     uint64_t genomeSizeRef) const {
    // Anchors must be the locations of the up- and downstream anchors, not the
    // data frame with ref and qry coordinates.
    // The scaling factor determines how fast the function falls when moving
    // away from an anchor.
    // Ideally, we define a half-life X_half, i.e. at a distance of X_half, the
    // model is at 0.5. With a scaling factor of 50 kb, X_half is at 20 kb (with
    // 100 kb at 10 kb).
    // score = 0.5^{minDist * genomeSizeRef / (genomeSize * halfLifeDistance_)}
    uint32_t const minDist(std::min(loc - upBound, downBound - loc));
    double const exp(((1.0d*minDist) / halfLifeDistance_)
                     * ((1.0d*genomeSizeRef) / genomeSize));
    double const score(std::pow(0.5d, exp));
    assert(0 <= score && score <= 1);
    return score;
}

void
Ipp::projectCoords(
    std::string const& refSpecies,
    std::string const& qrySpecies,
    std::vector<Coords> const& refCoords,
    unsigned const nThreads,
    OnProjectCoordsJobDoneCallback const& onJobDoneCallback) {
    // Calls projectCoord() on the given list of refCoords.
    // If nThreads > 1 then that many worker threads are started.
    // For each completed job the onJobDoneCallback() is called with the result.
    // The call to onJobDoneCallback() can come from any thread but no
    // concurrent calls will be made.

    std::mutex mutex;
    std::vector<Coords> jobs(refCoords.rbegin(), refCoords.rend());
    std::exception_ptr workerException;

    cancel_ = false;

    auto const worker = [&]() {
        while (true) {
            // Get the next job.
            Coords refCoord;
            {
                std::lock_guard const lockGuard(mutex);
                if (jobs.empty()) {
                    // All jobs done.
                    return;
                } else if (workerException) {
                    // An exception occured on a different thread. Abort.
                    return;
                } else if (cancel_) {
                    return;
                }

                refCoord = jobs.back();
                jobs.pop_back();
            }

            try {
                // Execute the next job (while not holding the mutex!).
                Ipp::CoordProjection const coordProjection(
                    projectCoord(refSpecies, qrySpecies, refCoord));

                // Call the callback (while holding the mutex).
                {
                    std::lock_guard const lockGuard(mutex);
                    onJobDoneCallback(refCoord, coordProjection);
                }
            } catch (...) {
                std::lock_guard const lockGuard(mutex);
                workerException = std::current_exception();
                return;
            }
        }
    };

    if (nThreads <= 1) {
        // Just execute the worker in this thread.
        worker();
    } else {
        // Create the threads.
        std::vector<std::thread> threads;
        for (unsigned i(0); i < nThreads; ++i) {
            threads.emplace_back(worker);
        }

        // Wait for the threads to complete.
        for (auto& thread : threads) {
            thread.join();
        }
    }

    // Forward any exception that occured in a worker.
    if (workerException) {
        std::rethrow_exception(workerException);
    }
}

void
Ipp::cancel() {
    // Cancel ongoing project_coords() call.
    cancel_ = true;
}

namespace {

struct ShortestPathMapKey {
    std::string species;
    Ipp::Coords coords;

    ShortestPathMapKey() {}
    ShortestPathMapKey(std::string const& species, Ipp::Coords const& coords)
        : species(species)
        , coords(coords)
    {}

    bool operator==(ShortestPathMapKey const& other) const {
        return species == other.species && coords == other.coords;
    }
};

struct ShortestPathMapEntry {
    double score;
    Ipp::Anchors anchors;
    ShortestPathMapKey const* prevKey;
    ShortestPathMapEntry const* prevEntry;

    ShortestPathMapEntry()
        : score(0) , prevKey(nullptr), prevEntry(nullptr)
    {}
    ShortestPathMapEntry(double score, 
                         Ipp::Anchors const& anchors,
                         ShortestPathMapKey const* prevKey,
                         ShortestPathMapEntry const* prevEntry)
        : score(score)
        , anchors(anchors)
        , prevKey(prevKey)
        , prevEntry(prevEntry)
    {}
};

struct OrangeEntry {
    double score;
    int pathLength;
    ShortestPathMapKey const* spKey;

    OrangeEntry(double score,
                int pathLength,
                ShortestPathMapKey const* spKey)
        : score(score)
        , pathLength(pathLength)
        , spKey(spKey)
    {}

    bool operator<(OrangeEntry const& other) const {
        // Use negative path length as shorter paths are better than longer
        // ones.
        return std::forward_as_tuple(score, -pathLength)
            < std::forward_as_tuple(other.score, -other.pathLength);
    }
};

} // namespace

template<>
struct std::hash<Ipp::Coords> {
    std::size_t operator()(Ipp::Coords const& coords) const noexcept {
        return std::hash<decltype(coords.chrom)>{}(coords.chrom)
            ^ (std::hash<decltype(coords.loc)>{}(coords.loc) << 1);
    }
};
template<>
struct std::hash<ShortestPathMapKey> {
    std::size_t operator()(ShortestPathMapKey const& spk) const noexcept {
        return std::hash<decltype(spk.species)>{}(spk.species)
            ^ (std::hash<decltype(spk.coords)>{}(spk.coords) << 1);
    }
};

Ipp::CoordProjection
Ipp::projectCoord(std::string const& refSpecies,
                  std::string const& qrySpecies,
                  Coords const& refCoords) const {
    bool const debug(false);
    if (debug) {
        std::cout.precision(16);
        std::cout << std::endl;
        std::cout << refSpecies << " " << qrySpecies << " "
                  << refCoords.chrom << ":" << refCoords.loc << std::endl;
    }

    uint64_t const genomeSizeRef(genomeSizes_.at(refSpecies));

    CoordProjection coordProjection;
    std::unordered_map<ShortestPathMapKey, ShortestPathMapEntry> shortestPath;

    ShortestPathMapKey const* const refSpKey(
        &shortestPath.try_emplace(ShortestPathMapKey(refSpecies, refCoords),
                                  1.0, Ipp::Anchors(), nullptr, nullptr)
        .first->first);

    std::priority_queue<OrangeEntry> orange; // greatest first.
    orange.emplace(1.0, 0, refSpKey);

    ShortestPathMapKey const* bestQrySpKey(nullptr);
    while (!orange.empty()) {
        OrangeEntry const current(orange.top());
        orange.pop();

        auto const it(shortestPath.find(*current.spKey));
        if (it != shortestPath.end() && it->second.score > current.score) {
            continue;
            // The current <species,coord> was already reached by a faster path,
            // ignore this path and go to the next species.
        }

        std::string const& currentSpecies(current.spKey->species);
        Coords const& currentCoords(current.spKey->coords);
        if (debug) {
            std::cout << "- " << currentSpecies << " " << current.score << " "
                      << currentCoords.chrom << ":" << currentCoords.loc
                      << std::endl;
        }
        
        if (currentSpecies == qrySpecies) {
            bestQrySpKey = current.spKey;
            break;
            // qry species reached, stop.
        }

        // Collect all the species that are traversed on the current path and
        // avoid visiting them again.
        std::unordered_set<std::string> speciesOnPath;
        ShortestPathMapKey const* currentSpKey(current.spKey);
        ShortestPathMapEntry const* currentSpEntry(&it->second);
        while (currentSpKey) {
            speciesOnPath.insert(currentSpKey->species);
            currentSpKey = currentSpEntry->prevKey;
            currentSpEntry = currentSpEntry->prevEntry;
        }

        for (auto const& nxt : pwalns_.at(currentSpecies)) {
            std::string const& nxtSpecies(nxt.first);

            // Don't visit a species twice on the same path.
            if (speciesOnPath.find(nxtSpecies) != speciesOnPath.end()) {
                continue;
            }

            if (debug) {
                std::cout << "--> " << nxtSpecies << std::endl;
            }

            std::optional<GenomicProjectionResult> const proj(
                projectGenomicLocation(currentSpecies,
                                       nxtSpecies,
                                       currentCoords,
                                       genomeSizeRef));
            if (!proj) {
                continue;
                // No path was found.
            }

            if(currentSpecies == refSpecies && nxtSpecies == qrySpecies) {
                // Direct projection.
                coordProjection.direct = *proj;
            }

            double const nxtScore(current.score * proj->score);
            const auto [it2, success] = shortestPath.try_emplace(
                {nxtSpecies, proj->nextCoords},
                nxtScore, proj->anchors, current.spKey, &it->second);
            if (success || it2->second.score < nxtScore) {
                if (!success) {
                    // There was already an entry in shortestPath but it had a
                    // worse score -> replace.
                    it2->second.score = nxtScore;
                    it2->second.anchors = proj->anchors;
                    it2->second.prevKey = current.spKey;
                    it2->second.prevEntry = &it->second;
                }

                // Only increase the path length if we don't reach the query
                // species as the next hop. This ensures that for the same
                // score we prefer the path that reaches the qry species first.
                int const nxtPathLength(nxtSpecies != qrySpecies
                                        ? current.pathLength + 1
                                        : current.pathLength);
                orange.emplace(nxtScore, nxtPathLength, &it2->first);
            }
        }
    }

    if (bestQrySpKey) {
        // Backtrace the shortest path from the reference to the given target
        // species (in reversed order).
        ShortestPathMapKey const* currentSpKey(bestQrySpKey);
        ShortestPathMapEntry const* currentSpEntry(
            &shortestPath.at(*currentSpKey));
        while (currentSpKey) {
            coordProjection.multiShortestPath.emplace_back(
                currentSpKey->species,
                currentSpKey->coords,
                currentSpEntry->score,
                currentSpEntry->anchors);
            currentSpKey = currentSpEntry->prevKey;
            currentSpEntry = currentSpEntry->prevEntry;
        }

        // Reverse the shortest path list to have it in the right order.
        std::reverse(coordProjection.multiShortestPath.begin(),
                     coordProjection.multiShortestPath.end());
    }

    return coordProjection;
}

namespace {

template<typename T>
void
updateChromCount(std::unordered_map<Ipp::ChromId, unsigned>& chromCount,
                 T const& v) {
    for (Ipp::PwalnEntry const* entry : v) {
        ++chromCount[entry->qryChrom()];
    }
}

template<typename T, typename... Args>
void
updateChromCount(std::unordered_map<Ipp::ChromId, unsigned>& chromCount,
                 T const& v,
                 Args const&... args) {
    updateChromCount(chromCount, v);
    updateChromCount(chromCount, args...);
}

template<typename... Args>
Ipp::ChromId
computeMajorChrom(Args const&... args) {
    std::unordered_map<Ipp::ChromId, unsigned> chromCount;
    updateChromCount(chromCount, args...);

    unsigned maxCount(0);
    Ipp::ChromId maxChrom(0);
    for (auto const& [chrom, count] : chromCount) {
        if (count > maxCount) {
            maxChrom = chrom;
            maxCount = count;
        }
    }
    return maxChrom;
}

template<typename T>
unsigned
insertIfMajorQryChromosome(
    T& anchors,
    Ipp::ChromId majorChrom,
    std::vector<Ipp::PwalnEntry const*>* closestAnchors) {
    unsigned numInserted(0);
    for (Ipp::PwalnEntry const* e : anchors) {
        if (e->qryChrom() == majorChrom) {
            closestAnchors->push_back(e);
            ++numInserted;
        }
    }
    return numInserted;
}

} // namespace

std::optional<Ipp::GenomicProjectionResult>
Ipp::projectGenomicLocation(std::string const& refSpecies,
                            std::string const& qrySpecies,
                            Coords const& refCoords,
                            uint64_t genomeSizeRef) const {
    auto const it1(pwalns_.find(refSpecies));
    if (it1 == pwalns_.end()) {
        // There is no pairwise alignment for the ref species.
        return {};
    }
    auto const it2(it1->second.find(qrySpecies));
    if (it2 == it1->second.end()) {
        // There is no pairwise alignment for the qry species.
        return {};
    }

    auto const anchors(
        getAnchors(it2->second, refSpecies, refCoords, qrySpecies));
    if (!anchors) {
        // If no or only one anchor is found because of border region, return 0
        // score and empty coordinate string.
        return {};
    }

    uint32_t const refLoc(refCoords.loc);

    // Compute the qryLoc by linear interpolation: Consider where refLoc lies
    // between the ref coords of the up- and downstream anchors and project that
    // to the qry coords of the anchors.
    // Note: The qry coords might be reversed (up.start > up.end). Either both
    //       anchors are reversed or both are not.
    //       The ref coords are never reversed.
    uint32_t qryLoc;
    double score;
    if (anchors->upstream == anchors->downstream) {
        // refLoc lies on an aligment.
        //  [  up.ref  ]
        //  [ down.ref ]
        //          x
        uint32_t const refUpBound(anchors->upstream.refStart());
        uint32_t const refDownBound(anchors->upstream.refEnd());
        uint32_t const qryUpBound(anchors->upstream.qryStart());

        assert(refUpBound <= refLoc && refLoc <= refDownBound);

        // Both endpoints are inclusive and the ref region is guaranteed to be
        // the same size than the qry region.
        score = 1.0d;
        qryLoc = !anchors->upstream.isQryReversed()
            ? qryUpBound + (refLoc - refUpBound)
            : qryUpBound - (refLoc - refUpBound);
    } else {
        // [ up.ref ]  x    [ down.ref ]
        assert(anchors->upstream.isQryReversed()
               == anchors->downstream.isQryReversed());
        uint32_t const refUpBound(anchors->upstream.refEnd());
        uint32_t const refDownBound(anchors->downstream.refStart()); 
        uint32_t const qryUpBound(anchors->upstream.qryEnd());
        uint32_t const qryDownBound(anchors->downstream.qryStart());

        // Both endpoints are exclusive.
        assert(refUpBound < refLoc && refLoc < refDownBound);

        // ONLY USE DISTANCE TO CLOSE ANCHOR AT REF SPECIES, because at the qry
        // species it should be roughly the same as it is a projection of the
        // reference.
        score = projectionScore(refLoc,
                                refUpBound,
                                refDownBound,
                                genomeSizes_.at(refSpecies),
                                genomeSizeRef);

        // +0.5 to bring the projection into the middle of the projected qry
        // region of potentially different size:
        //     up.refEnd = 10, down.refStart = 20
        //     up.qryEnd = 110, down.qryStart = 150
        //     refLoc = 11
        //     qryLoc = 110 + ((11 - 10 + 0.5) / (20-10)) * (150-110)
        //            = 110 + 1.5/10 * 40 = 116
        //     (vs. 114 w/o the +0.5).
        double const relativeRefLoc(
            1.0d*(refLoc - refUpBound + 0.5) / (refDownBound - refUpBound));
        bool const isQryReversed(anchors->upstream.isQryReversed());
        qryLoc = !isQryReversed
            ? qryUpBound + relativeRefLoc*(qryDownBound - qryUpBound)
            : qryUpBound - relativeRefLoc*(qryUpBound - qryDownBound);

        // <= and >= below in case the qry range is of size <= 1 (then the
        // projection is onto the lower boundary).
        assert((!isQryReversed && qryUpBound <= qryLoc && qryLoc < qryDownBound)
               ||(isQryReversed && qryUpBound > qryLoc && qryLoc >= qryDownBound));
    }

    return {{score, {anchors->upstream.qryChrom(), qryLoc}, *anchors}};
}

std::optional<Ipp::Anchors>
Ipp::getAnchors(Pwaln const& pwaln,
                std::string const& refSpecies,
                Coords const& refCoords,
                std::string const& qrySpecies) const {
    // First define anchors upstream, downstream and ovAln, then do major-chrom
    // and collinearity test, then either return overlapping anchor or closest
    // anchors.
    // Take orientation into account for the anchor definition. If start > end,
    // then the aln is to the '-' strand.
    // For speed reasons only select the first `topn` entries. The rest just
    // takes longer to compute min / max and most likely will (and should) not
    // be an anchor anyways.

    // Test collinearity of anchors: take top 20 in each direction (top 10 
    // produced many locally collinear pwalns that were still non-collinear
    // outliers in the global view of the GRB).
    // Note: using ungapped chain blocks might require n to be even larger.
    unsigned const minn(5);
    unsigned const topn(20);

    uint32_t const refLoc(refCoords.loc);

    auto const pwalnEntriesIt(pwaln.find(refCoords.chrom));
    if (pwalnEntriesIt == pwaln.end()) {
        // No pwaln entry for this refCoords.chrom.
        return {};
    }
    std::vector<PwalnEntry> const& pwalnEntries(pwalnEntriesIt->second);

    // Find the topn entries by largest(smallest) refEnd(refStart) in the
    // upstream(downstream) anchors.

    // Binary search for the closest upstream anchor (the first with
    // refStart > refLoc).
    auto const closestDownstreamAnchorIt(std::upper_bound(
            pwalnEntries.begin(),
            pwalnEntries.end(),
            refLoc,
            [](uint32_t loc, PwalnEntry const& e) {
                return loc < e.refStart();
            }));

    // Find the downstream anchors.
    std::vector<PwalnEntry const*> anchorsDownstream;
    anchorsDownstream.reserve(topn);
    for (auto it(closestDownstreamAnchorIt); it != pwalnEntries.end(); ++it) {
        // downstream anchor
        //    x     [ anchor ]
        PwalnEntry const& pwalnEntry(*it);
        assert(refLoc < pwalnEntry.refStart());
        anchorsDownstream.push_back(&pwalnEntry);
        if(anchorsDownstream.size() == topn) {
            // We found the topn closest anchors with refStart > refLoc.
            // Since the pwalnEntries are sorted by refStart, all the
            // anchors to come are further away than what we have already
            // seen.
            break;
        }
    }

    // Find the ovAln and upstream anchors. Walk backwards on the pwalnEntries
    // starting from closestDownstreamAnchorIt and walk until
    // e.refStart+maxAnchorLength_ < furthestUpstreamAnchor.refEnd.
    // For this, keep a priority queue of the furthest upstream anchor.
    auto const compGreaterRefEnd = [](PwalnEntry const* lhs,
                                      PwalnEntry const* rhs) {
        return std::forward_as_tuple(lhs->refEnd(),
                                     lhs->refStart(),
                                     lhs->qryStart(),
                                     lhs->qryChrom())
            > std::forward_as_tuple(rhs->refEnd(),
                                    rhs->refStart(),
                                    rhs->qryStart(),
                                    rhs->qryChrom());
    };
    std::vector<PwalnEntry const*> anchorsUpstreamPqContainer;
    anchorsUpstreamPqContainer.reserve(topn+1);
    std::priority_queue<PwalnEntry const*,
                        decltype(anchorsUpstreamPqContainer),
                        decltype(compGreaterRefEnd)> anchorsUpstreamPq(
            compGreaterRefEnd,
            std::move(anchorsUpstreamPqContainer));
    std::vector<PwalnEntry const*> ovAln;
    for (auto it(closestDownstreamAnchorIt); it-- != pwalnEntries.begin();) {
        // Walk upstream on the chromosome.
        PwalnEntry const& pwalnEntry(*it);
        if (pwalnEntry.refEnd() < refLoc) { // refEnd is inclusive
            // upstream anchor
            // [ anchor ]    x
            if (anchorsUpstreamPq.size() == topn) {
                if (pwalnEntry.refStart() + maxAnchorLength_
                    < anchorsUpstreamPq.top()->refEnd()) {
                    // This anchor is so far away from the currently furthest
                    // upstream anchor that there is no further pwalnEntry2
                    // (with pwalnEntry2.refStart < pwalnEntry.refStart) that
                    // has pwalnEntry2.refEnd > pwalnEntry.refEnd.
                    break;
                }

                if (pwalnEntry.refEnd() < anchorsUpstreamPq.top()->refEnd()) {
                    // Prevent adding an entry that would immediately be removed
                    // again.
                    continue;
                }
            }

            anchorsUpstreamPq.push(&pwalnEntry);
            if (anchorsUpstreamPq.size() > topn) {
                // Remove surplus anchors that are too far away.
                anchorsUpstreamPq.pop();
            }
        } else {
            // refLoc lies on an alignment block.
            // [ anchor ]
            //      x
            assert(pwalnEntry.refStart() <= refLoc
                   && refLoc <= pwalnEntry.refEnd());
            ovAln.push_back(&pwalnEntry);
        }
    }
    assert(anchorsUpstreamPq.size() <= topn);

    // Convert the anchorsUpstream priority queue into a vector.
    std::vector<PwalnEntry const*> anchorsUpstream;
    anchorsUpstream.reserve(anchorsUpstreamPq.size());
    while (!anchorsUpstreamPq.empty()) {
        anchorsUpstream.push_back(anchorsUpstreamPq.top());
        anchorsUpstreamPq.pop();
    }

    // MAJOR CHROMOSOME: Retain anchors that point to the majority chromosome in
    // top n of both up- and downstream anchors.
    ChromId const majorChrom(
        computeMajorChrom(ovAln, anchorsUpstream, anchorsDownstream));
    std::vector<PwalnEntry const*> closestAnchors;
    closestAnchors.reserve(
        anchorsUpstream.size() + ovAln.size() + anchorsDownstream.size());
    unsigned const numUpstream(insertIfMajorQryChromosome(anchorsUpstream,
                                                          majorChrom,
                                                          &closestAnchors));
    insertIfMajorQryChromosome(ovAln, majorChrom, &closestAnchors);
    unsigned const numDownstream(insertIfMajorQryChromosome(anchorsDownstream,
                                                            majorChrom,
                                                            &closestAnchors));
    if (!numUpstream || !numDownstream) {
        // Require minimum of 1 anchor on each side. Later, the minimum total
        // number of collinear anchors will be set to `minn` (but one side is
        // allowed to have as little as 1 anchor).
        return {};
    }

    // Sort the closestAnchors entries by increasing refStart. That is necessary
    // as the anchorsUpstream were previously sorted by decreasing refEnd.
    auto const compLessRefStart = [](PwalnEntry const* lhs,
                                     PwalnEntry const* rhs) {
        return std::forward_as_tuple(lhs->refStart(),
                                     lhs->refEnd(),
                                     lhs->qryStart(),
                                     lhs->qryChrom())
            < std::forward_as_tuple(rhs->refStart(),
                                    rhs->refEnd(),
                                    rhs->qryStart(),
                                    rhs->qryChrom());
    };
    std::sort(closestAnchors.begin(), closestAnchors.end(), compLessRefStart);

    // COLLINEARITY: Remove pwalns pointing to outliers by getting the longest
    // sorted subsequence of the top n of both up- and downstream anchors.
    // Compute longest collinear anchor subsequence (while considering both
    // normal and reversed direction).
    closestAnchors = longestSubsequence(closestAnchors);

    // Set minimum number of collinear anchors to `minn` (for species pairs with
    // very large evol. distances setting a lower boundary for the number of
    // collinear anchors will help reduce false positives).
    if (closestAnchors.size() < minn) {
        return {};
    }
  
    // Check if the original ovAln is still present (or ever was) in the
    // filtered closestAnchors (that include a potential ovAln);
    // if not, it was an outlier alignment and was filtered out
    PwalnEntry const* closestUpstreamAnchor(nullptr);
    PwalnEntry const* closestOvAlnAnchor(nullptr);
    PwalnEntry const* closestDownstreamAnchor(nullptr);
    for (PwalnEntry const* anchor : closestAnchors) {
        if (anchor->refEnd() < refLoc) {
            if (!closestUpstreamAnchor
                || closestUpstreamAnchor->refEnd() < anchor->refEnd()) {
                closestUpstreamAnchor = anchor;
            }
        } else if (refLoc < anchor->refStart()) {
            if (!closestDownstreamAnchor
                || anchor->refStart() < closestDownstreamAnchor->refStart()) {
                closestDownstreamAnchor = anchor;

                // The anchors that follow this one will only be worse.
                break;
            }
        } else {
            if (false&&closestOvAlnAnchor) {
                // An overlapping direct mapping found.
                std::cerr << std::endl;
                std::cerr << format("WARNING: Overlapping direct mapping for "
                                    "(%s -> %s @ %s:%u): ",
                                    refSpecies.c_str(),
                                    qrySpecies.c_str(),
                                    chroms_.at(refCoords.chrom).c_str(),
                                    refCoords.loc) << std::endl;

                auto const printAnchor = [this](PwalnEntry const* anchor) {
                    std::cerr << "    refStart: " << anchor->refStart() << std::endl;
                    std::cerr << "    refEnd: " << anchor->refEnd() << std::endl;
                    std::cerr << "    qryChrom: " << chroms_.at(anchor->qryChrom()) << std::endl;
                    std::cerr << "    qryStart: " << anchor->qryStart() << std::endl;
                    std::cerr << "    qryEnd: " << anchor->qryEnd() << std::endl;
                };
                std::cerr << "  anchor 1" << std::endl;
                printAnchor(closestOvAlnAnchor);
                std::cerr << "  anchor 2" << std::endl;
                printAnchor(anchor);
                std::cerr << std::endl;
            }
            closestOvAlnAnchor = anchor;
        }
    }
    if (closestOvAlnAnchor) {
        // We found a direct mapping.
        return {{*closestOvAlnAnchor, *closestOvAlnAnchor}};
    } else {
        if (!closestUpstreamAnchor || !closestDownstreamAnchor) {
            // Not both up- and downstream anchors were found (e.g. at synteny
            // break points where one side does not have any anchors to the
            // majority chromosome)
            return {};
        }

        return {{*closestUpstreamAnchor, *closestDownstreamAnchor}};
    }
}

namespace {

std::vector<Ipp::PwalnEntry const*>
longestSubsequence(std::vector<Ipp::PwalnEntry const*> const& seq,
                   std::function<bool(Ipp::PwalnEntry const*)> const& filter,
                   std::function<int(Ipp::PwalnEntry const*)> const& qryStart,
                   std::function<int(Ipp::PwalnEntry const*)> const& qryEnd) {
    // Finds the longest strictly increasing subsequence (in regards to the
    // given greater-than function). Only elements for which filter(seq[i])
    // returns true are considered.
    // O(n log k) algorithm.

    // m[i] contains the index to the smallest value in seq[] that is the end of
    // a subsequence of length i+1.
    std::vector<unsigned> m;
    m.reserve(seq.size());

    if (!seq.size()) {
        return {};
    }

    // prev[i] contains the index of the element in seq that is the one before
    // seq[i] in the longest subsequence for seq[i].
    std::vector<unsigned> prev(seq.size());

    for (unsigned i(0); i < seq.size(); ++i) {
        // If the next element seq[i] is greater than the last element of the
        // current longest subsequence seq[m.back()], just push it to the end of
        // `m` and continue.
        if (!filter(seq[i])) {
            continue;
        }

        if (!m.size()) {
            // This is the first element that matches the filter. Just add it to
            // m.
            m.push_back(i);
            continue;
        }

        if (qryEnd(seq[m.back()]) < qryStart(seq[i])) {
            prev[i] = m.back();
            m.push_back(i);
            continue;
        }

        // Binary search to find the smallest element referenced by m which is
        // just bigger than seq[i].
        // Note : Binary search is performed on m (and not seq).
        // Size of m is always <=i and hence contributes O(log i) to the
        // complexity.
        unsigned u(0);
        unsigned v(m.size()-1);
        while(u < v) {
            unsigned const mid((u + v) / 2);
            if (qryEnd(seq[m[mid]]) < qryStart(seq[i])) {
                u = mid+1;
            } else {
                v = mid;
            }
        }

        // Update m if the new value is smaller than the previously referenced
        // one.
        if (qryEnd(seq[i]) < qryEnd(seq[m[u]])) {
            if (u > 0) {
                prev[i] = m[u-1];
            }
            m[u] = i;
        }
    }

    // Backtrace the longest subsequence into res.
    std::vector<Ipp::PwalnEntry const*> res(m.size());
    unsigned v(m.back());
    for (unsigned u(m.size()); u; --u) {
        res[u-1] = seq[v];
        v = prev[v];
    }
    return res;
}

} // namespace

std::vector<Ipp::PwalnEntry const*>
Ipp::longestSubsequence(std::vector<PwalnEntry const*> const& seq) {
    std::vector<PwalnEntry const*> const inc(::longestSubsequence(
            seq,
            [](PwalnEntry const* e) {
                return !e->isQryReversed();
            },
            [](PwalnEntry const* e) {
                return e->qryStart();
            },
            [](PwalnEntry const* e) {
                return e->qryEnd();
            }));

    std::vector<PwalnEntry const*> const dec(::longestSubsequence(
            seq,
            [](PwalnEntry const* e) {
                return e->isQryReversed();
            },
            [](PwalnEntry const* e) {
                return -1 * (int)e->qryStart();
            },
            [](PwalnEntry const* e) {
                return -1 * (int)e->qryEnd();
            }));

    // Sanity check: The entries in the inc/dec list should be strictly
    // increasing/decreasing.
    uint32_t loc(0);
    for (PwalnEntry const* e : inc) {
        assert((!loc && !e->qryStart()) || loc < e->qryStart());
        assert(e->qryStart() <= e->qryEnd());
        loc = e->qryEnd();
    }
    loc = std::numeric_limits<uint32_t>::max();
    for (PwalnEntry const* e : dec) {
        assert(loc > e->qryEnd());
        assert(e->qryStart() >= e->qryEnd());
        loc = e->qryStart();
    }

    return inc.size() >= dec.size() ? inc : dec;
}
