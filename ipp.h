#pragma once

#include <functional>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class Ipp {
public:
    struct PwalnEntry {
        // Attention: Keep the below order of the members for alignment reasons!
        uint32_t refStart;
        uint32_t refEnd;
        uint32_t qryStart;
        uint32_t qryEnd;
        uint16_t refChrom;
        uint16_t qryChrom;

        PwalnEntry() {}
        PwalnEntry(PwalnEntry const& other)
            : refStart(other.refStart)
            , refEnd(other.refEnd)
            , qryStart(other.qryStart)
            , qryEnd(other.qryEnd)
            , refChrom(other.refChrom)
            , qryChrom(other.qryChrom)
        {}

        bool isQryReversed() const {
            return qryStart > qryEnd;
        }

        bool operator==(PwalnEntry const& other) const {
            return refChrom == other.refChrom
                && refStart == other.refStart
                && refEnd == other.refEnd
                && qryChrom == other.qryChrom
                && qryStart == other.qryStart
                && qryEnd == other.qryEnd;
        }
    };
    using Pwaln = std::unordered_map<uint16_t, std::vector<PwalnEntry>>;
    using Pwalns = std::unordered_map<std::string, std::unordered_map<std::string, Pwaln>>;
    // The map of pairwise alignments: [sp1][sp2][ref_chrom] -> [PwalnEntry]
    // The entries in the vector are sorted by [refStart, qryChrom, qryStart].
    // refChrom is the key and also in the entry as PwalnEntry has enough space
    // to keep it (due to alignment) and that way we can bulk-read the data in
    // loadPwalns().

    void loadPwalns(std::string const& fileName);
    // Reads the chromosomes and pwalns from the given file.

    void loadGenomeSizes(std::string const& dirName);
    // Reads the genome sizes from the files in the given directory.

    void setHalfLifeDistance(unsigned halfLifeDistance);
    // Sets the half-life distance.

    uint16_t chromIdFromName(std::string const& chromName) const;
    // Looks up the given chromosome name in chroms_ and returns its id.

    std::string const& chromName(uint16_t chromId) const;
    // Returns the name of the chromosome with the given id.

    struct Coords {
        uint16_t chrom;
        uint32_t loc;

        Coords() {}
        Coords(uint16_t chrom, uint32_t loc) : chrom(chrom), loc(loc) {}

        bool operator<(Coords const& other) const {
            return std::tie(chrom, loc) < std::tie(other.chrom, other.loc);
        }
    };

    struct Anchors {
        PwalnEntry upstream;
        PwalnEntry downstream;

        Anchors() {}
        Anchors(PwalnEntry upstream, PwalnEntry downstream)
            : upstream(upstream)
            , downstream(downstream)
        {}
    };

    struct ShortestPathEntry {
        double score;
        std::string prevSpecies;
        Coords coords;
        Anchors anchors;

        ShortestPathEntry() {}
        ShortestPathEntry(double score, 
                          std::string const& prevSpecies,
                          Coords const& coords,
                          Anchors const& anchors)
            : score(score)
            , prevSpecies(prevSpecies)
            , coords(coords)
            , anchors(anchors)
        {}
    };
    using ShortestPath = std::unordered_map<std::string, ShortestPathEntry>;

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
        unsigned const nCores,
        OnProjectCoordsJobDoneCallback const& onJobDoneCallback) const;

private:
    CoordProjection projectCoord(std::string const& refSpecies,
                                 std::string const& qrySpecies,
                                 Coords const& refCoords) const;

    std::optional<GenomicProjectionResult> projectGenomicLocation(
        std::string const& refSpecies,
        std::string const& qrySpecies,
        Coords const& refCoords,
        double scalingFactor) const;

    std::optional<Anchors> getAnchors(Pwaln const& pwaln,
                                      Coords const& refCoords) const;

    static std::vector<PwalnEntry> longestSubsequence(
        std::vector<PwalnEntry> const& seq);
    // Searches the longest strictly increasing or decreasing subsequence of seq
    // and puts it in res.
    // O(n log k) algorithm.

    double getScalingFactor(unsigned genomeSize) const;

    double projectionScore(uint32_t loc,
                           uint32_t leftBound,
                           uint32_t rightBound,
                           unsigned genomeSize,
                           double scalingFactor) const;

private:
    std::vector<std::string> chroms_;
    Pwalns pwalns_;
    std::unordered_map<std::string, unsigned> genomeSizes_;
    unsigned halfLifeDistance_;
};


