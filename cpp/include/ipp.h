#pragma once

#include <functional>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>

class Ipp {
public:
    using ChromId = uint32_t;

    class PwalnEntry {
    public:
        PwalnEntry()
            : refStart_(0)
            , qryStart_(0)
            , qryChrom_(0)
            , lengthAndStrand_(0)
        {}
        PwalnEntry(PwalnEntry const& other)
            : refStart_(other.refStart_)
            , qryStart_(other.qryStart_)
            , qryChrom_(other.qryChrom_)
            , lengthAndStrand_(other.lengthAndStrand_)
        {}

        uint32_t refStart() const {
            return refStart_;
        }
        uint32_t refEnd() const {
            // Inclusive. Ref is always increasing no matter the strand.
            return refStart_ + length() - 1;
        }

        ChromId qryChrom() const {
            return qryChrom_;
        }
        uint32_t qryStart() const {
            return qryStart_;
        }
        uint32_t qryEnd() const {
            return !isQryReversed()
                ? qryStart_ + length() - 1
                : qryStart_ - length() + 1;
        }

        uint16_t length() const {
            // Remove the MSB bit from lengthAndStrand.
            return lengthAndStrand_ & ~(1<<15);
        }
        bool isQryReversed() const {
            // Returns true if the MSB of the lengthAndStrand is 1.
            return lengthAndStrand_ & (1<<15);
        }

        bool operator==(PwalnEntry const& other) const {
            return refStart_ == other.refStart_
                && qryStart_ == other.qryStart_
                && qryChrom_ == other.qryChrom_
                && lengthAndStrand_ == other.lengthAndStrand_;
        }

    private:
        // Attention: Keep the below order of the members for alignment reasons!
        uint32_t refStart_;
        uint32_t qryStart_;
        ChromId qryChrom_;

        // MSB is 1 for negative strand.
        // Note: This is not the two's complement so don't treat this as a
        // signed integer.
        uint16_t lengthAndStrand_;
    };
    using Pwaln = std::unordered_map<ChromId, std::vector<PwalnEntry>>;
    using Pwalns = std::unordered_map<std::string, std::unordered_map<std::string, Pwaln>>;
    // The map of pairwise alignments: [sp1][sp2][ref_chrom] -> [PwalnEntry]
    // The entries in the vector are sorted by [refStart, qryChrom, qryStart].
    // refChrom is the key and also in the entry as PwalnEntry has enough space
    // to keep it (due to alignment) and that way we can bulk-read the data in
    // loadPwalns().

    void loadPwalns(std::string const& fileName);
    // Reads the chromosomes and pwalns from the given file.

    uint64_t getGenomeSize(std::string const& speciesName);
    // Returns the genome size for a given species name

    void setHalfLifeDistance(unsigned halfLifeDistance);
    // Sets the half-life distance.

    std::optional<ChromId> chromIdFromName(std::string const& chromName) const;
    // Looks up the given chromosome name in chroms_ and returns its id.

    std::string const& chromName(ChromId chromId) const;
    // Returns the name of the chromosome with the given id.

    struct Coords {
        ChromId chrom;
        uint32_t loc;

        Coords()
            : chrom(0)
            , loc(0)
        {}
        Coords(ChromId chrom, uint32_t loc) : chrom(chrom), loc(loc) {}

        bool operator<(Coords const& other) const {
            return std::tie(chrom, loc) < std::tie(other.chrom, other.loc);
        }
        bool operator==(Coords const& other) const {
            return chrom == other.chrom
                && loc == other.loc;
        }
    };

    struct Anchors {
        PwalnEntry upstream;
        PwalnEntry downstream;

        Anchors() {}
        Anchors(PwalnEntry const& upstream, PwalnEntry const& downstream)
            : upstream(upstream)
            , downstream(downstream)
        {}
    };

    struct ShortestPathEntry {
        std::string species;
        Coords coords;
        double score;
        Anchors anchors;

        ShortestPathEntry(std::string const& species,
                          Coords const& coords,
                          double score,
                          Anchors const& anchors)
            : species(species)
            , coords(coords)
            , score(score)
            , anchors(anchors)
        {}
    };
    using ShortestPath = std::vector<ShortestPathEntry>;

    struct GenomicProjectionResult {
        double score;
        Coords nextCoords;
        Anchors anchors;

        GenomicProjectionResult() {}
        GenomicProjectionResult(double score,
                                Coords const& nextCoords,
                                Anchors const& anchors)
            : score(score)
            , nextCoords(nextCoords)
            , anchors(anchors)
        {}
    };

    struct CoordProjection {
        std::optional<GenomicProjectionResult> direct;
        ShortestPath multiShortestPath;
    };

    using OnProjectCoordsJobDoneCallback =
        std::function<void(Ipp::Coords const&, Ipp::CoordProjection const&)>;

    void projectCoords(
        std::string const& refSpecies,
        std::string const& qrySpecies,
        std::vector<Coords> const& refCoords,
        unsigned const nThreads,
        OnProjectCoordsJobDoneCallback const& onJobDoneCallback);
    // Calls projectCoord() on the given list of refCoords.
    // If nThreads > 1 then that many worker processes are started.
    // For each completed job the onJobDoneCallback() is called with the result.
    // The call to onJobDoneCallback() can come from any thread but no
    // concurrent calls will be made.

    void cancel();
    // Cancel ongoing project_coords() call.

private:
    CoordProjection projectCoord(std::string const& refSpecies,
                                 std::string const& qrySpecies,
                                 Coords const& refCoords) const;

    std::vector<GenomicProjectionResult> projectGenomicLocation(
        std::string const& refSpecies,
        std::string const& qrySpecies,
        Coords const& refCoords,
        uint64_t genomeSizeRef) const;

    std::vector<Anchors> getAnchors(Pwaln const& pwaln,
                                    std::string const& refSpecies,
                                    Coords const& refCoords,
                                    std::string const& qrySpecies) const;

    static std::vector<PwalnEntry const*> longestSubsequence(
        std::vector<PwalnEntry const*> const& seq);
    // Searches the longest strictly increasing or decreasing subsequence of seq
    // and puts it in res.
    // O(n log k) algorithm.

    double projectionScore(uint32_t loc,
                           uint32_t upBound,
                           uint32_t downBound,
                           uint64_t genomeSize,
                           uint64_t genomeSizeRef) const;

private:
    std::vector<std::string> chroms_;
    Pwalns pwalns_;
    std::unordered_map<std::string, uint64_t> genomeSizes_;
    unsigned halfLifeDistance_;
    uint16_t maxAnchorLength_;
    volatile bool cancel_;
};

template <typename ...Args>
std::string
format(char const* fmt, Args&& ...args) {
    auto const len(std::snprintf(nullptr, 0, fmt, std::forward<Args>(args)...));

    std::string ret(len+1, '\0');
    std::sprintf(ret.data(), fmt, std::forward<Args>(args)...);
    return ret;
}
