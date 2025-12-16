/******************************************************************************
 * Copyright (c) 2025, Guilhem Villemin (guilhem.villemin@altametris.com)
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include "MtaProcessingFilter.hpp"

#include <pdal/pdal_internal.hpp>
#include <pdal/util/ProgramArgs.hpp>
#include <pdal/SpatialReference.hpp>
#include <cmath>
#include <limits>
#include <algorithm>


namespace pdal
{



static StaticPluginInfo const s_info
{
    "filters.mta",
    "MTA (Multiple-Time-Around) zone analysis and correction filter",
    "http://pdal.io/stages/filters.mta.html"
};

CREATE_STATIC_STAGE(MtaProcessingFilter, s_info)

std::string MtaProcessingFilter::getName() const
{
    return s_info.name;
}


struct ShotData {
    std::vector<size_t> echoes;
};

MtaProcessingFilter::MtaProcessingFilter()
    : m_mode("auto")
    , m_fixedZone(1)
    , m_gateNear(0.0f)
    , m_gateFar(0.0f)
    , m_maxZoneCandidates(4)
    , m_unambiguousRangeOverride(0.0)
    , m_maxOriginDistance(0.001)
    , m_maxAngleDifference(0.1)
    , m_jumpTolerance(0.05)
{}

MtaProcessingFilter::~MtaProcessingFilter()
{}

void MtaProcessingFilter::addArgs(ProgramArgs& args)
{
    args.add("mode", "MTA correction mode: fixed, rangegate, auto", 
             m_mode, "auto");
    args.add("fixed_zone", "Fixed zone number for 'fixed' mode", m_fixedZone, uint8_t(1));
    args.add("gate_near", "Near gate boundary in meters for 'rangegate' mode", m_gateNear, 0.0f);
    args.add("gate_far", "Far gate boundary in meters for 'rangegate' mode", m_gateFar, 0.0f);
    args.add("max_zones", "Maximum number of zone candidates to evaluate", 
             m_maxZoneCandidates, uint8_t(4));
    args.add("unambiguous_range", 
             "Unambiguous range override in meters. If set, overrides the UnambiguousRange dimension. "
             "Otherwise, reads from dimension.",
             m_unambiguousRangeOverride, 0.0);
        args.add("max_origin_distance",
                 "Maximum distance in meters between beam origins for matching (moving platform tolerance)",
                 m_maxOriginDistance, 0.001);
        args.add("max_angle_difference",
                 "Maximum angle difference in degrees between beam directions for matching",
                 m_maxAngleDifference, 0.1);
    args.add("jump_tolerance",
             "Tolerance (meters) around multiples of unambiguous range for acceptable jumps between consecutive shots",
             m_jumpTolerance, 0.05);
}

void MtaProcessingFilter::initialize()
{
    // Validate mode
    if (m_mode != "fixed" && m_mode != "rangegate" && m_mode != "auto")
    {
        throwError("Invalid MTA mode '" + m_mode + "'. Must be one of: fixed, rangegate, auto");
    }
    
    // Validate parameters for specific modes
    if (m_mode == "fixed" && m_fixedZone == 0)
    {
        throwError("fixed_zone must be >= 1 for 'fixed' mode");
    }
    
    if (m_mode == "rangegate" && m_gateFar <= m_gateNear)
    {
        throwError("gate_far must be > gate_near for 'rangegate' mode");
    }
}

void MtaProcessingFilter::prepared(PointTableRef table)
{
    (void)table; // no retained state; using layout directly from views
    
    // Verify we have ECEF coordinates (required for all modes)
    // Even fixed/rangegate need consistent reference frame for position correction
    const SpatialReference& srs = table.spatialReference();
    
    if (!srs.empty())
    {
        // Check if SRS is ECEF (EPSG:4978)
        SpatialReference ecefRef("EPSG:4978");
        if (!srs.equals(ecefRef))
        {
            std::string name = srs.getName();
            std::string msg = std::string("MTA processing requires ECEF coordinates (EPSG:4978) for accurate spatial analysis.\n") +
                std::string("Current SRS: ") + name + std::string("\n") +
                std::string("Projected or geographic coordinates will produce incorrect results.\n") +
                std::string("Use filters.georeference with coordinate_system=EPSG:4978 before this filter.");
            throwError(msg);
        }
    }
    else
    {
        log()->get(LogLevel::Info) << "MTA mode '" << m_mode 
            << "': No SRS defined, assuming ECEF coordinates." << std::endl;
    }

    // Verify UnambiguousRange availability or override for modes that require it
    if (m_mode != "fixed")
    {
        if (m_unambiguousRangeOverride <= 0.0)
        {
            // When no override is set, require the UnambiguousRange dimension in the layout
            const auto* layout = table.layout();
            bool hasUA = layout && layout->hasDim(Dimension::Id::UnambiguousRange);
            if (!hasUA)
            {
                throwError("UnambiguousRange is required for MTA processing in '" + m_mode +
                           "' mode. Provide the dimension or set 'unambiguous_range' override.");
            }
        }
    }
}

bool MtaProcessingFilter::processOne(PointRef& point)
{
    // Only used for "fixed" and "rangegate" modes
    
    using namespace Dimension;
    
    double measuredRange = point.getFieldAs<double>(Id::EchoRange);
    double unambiguousRange = 0.0;
    double beamOriginX = point.getFieldAs<double>(Id::BeamOriginX);
    double beamOriginY = point.getFieldAs<double>(Id::BeamOriginY);
    double beamOriginZ = point.getFieldAs<double>(Id::BeamOriginZ);
    double beamDirectionX = point.getFieldAs<double>(Id::BeamDirectionX);
    double beamDirectionY = point.getFieldAs<double>(Id::BeamDirectionY);
    double beamDirectionZ = point.getFieldAs<double>(Id::BeamDirectionZ);
    
    uint8_t rawZone = 1;
    
    if (m_mode == "fixed")
    {
        uint8_t correctedZone = m_fixedZone;
        if (correctedZone != rawZone && unambiguousRange > 0)
        {
            double correction = (correctedZone - 1) * unambiguousRange;
            double correctedRange = measuredRange + correction;
            point.setField(Id::X, beamOriginX + beamDirectionX * correctedRange);
            point.setField(Id::Y, beamOriginY + beamDirectionY * correctedRange);
            point.setField(Id::Z, beamOriginZ + beamDirectionZ * correctedRange);
        }
        return true;
    }
    else if (m_mode == "rangegate")
    {
        uint8_t correctedZone = correctRangeGate(measuredRange, unambiguousRange);
        if (correctedZone > 0 && unambiguousRange > 0)
        {
            double correction = (correctedZone - 1) * unambiguousRange;
            double correctedRange = measuredRange + correction;
            point.setField(Id::X, beamOriginX + beamDirectionX * correctedRange);
            point.setField(Id::Y, beamOriginY + beamDirectionY * correctedRange);
            point.setField(Id::Z, beamOriginZ + beamDirectionZ * correctedRange);
        }
        return true;
    }
    
    // Should not reach here for multi-round modes
    return false;
}

PointViewSet MtaProcessingFilter::run(PointViewPtr view)
{
    using namespace Dimension;

    // For simple streaming modes, just iterate and call processOne
    if (m_mode == "fixed" || m_mode == "rangegate")
    {
        PointRef pt = view->point(0);
        for (PointId idx = 0; idx < view->size(); ++idx)
        {
            pt.setPointId(idx);
            processOne(pt);
        }
        PointViewSet viewSet;
        viewSet.insert(view);
        return viewSet;
    }

    PointRef pt = view->point(0);
    // Group indices by shot (one vector of indices per shot)
    std::vector<std::vector<size_t>> shotGroups;
    std::vector<double> shotTimestamps;
    double lastShotTime = std::numeric_limits<double>::quiet_NaN();
    std::vector<size_t> currentGroup;
    for (size_t idx = 0; idx < view->size(); ++idx)
    {
        pt.setPointId(idx);
        double shotTimestamp = pt.getFieldAs<double>(Id::GpsTime);
        if (shotTimestamps.empty() || shotTimestamp != lastShotTime)
        {
            if (!shotTimestamps.empty())
                shotGroups.push_back(currentGroup);
            shotTimestamps.push_back(shotTimestamp);
            currentGroup.clear();
            lastShotTime = shotTimestamp;
        }
        currentGroup.push_back(idx);
    }
    if (!currentGroup.empty() || shotGroups.empty())
        shotGroups.push_back(currentGroup);

    // Remove empty shots
    for (auto& group : shotGroups)
    {
        if (group.size() == 1)
        {
            pt.setPointId(group[0]);
            uint8_t nret = pt.getFieldAs<uint8_t>(Id::NumberOfReturns);
            if (nret == 0)
                group.clear();
        }
    }

    if (shotGroups.empty())
    {
        PointViewSet viewSet;
        return viewSet;
    }

    // Windowed Viterbi zone assignment
    const size_t windowSize = 2 * m_maxZoneCandidates;
    std::vector<uint8_t> zones(shotGroups.size(), 1);
    for (size_t start = 0; start < shotGroups.size(); start += windowSize)
    {
        size_t end = std::min(start + windowSize, shotGroups.size());
        size_t winLen = end - start;
        if (winLen == 0) continue;

        // Build window of ShotData
        std::vector<ShotData> window;
        for (size_t i = start; i < end; ++i)
        {
            ShotData shot;
            shot.echoes = shotGroups[i];
            window.push_back(shot);
        }

        // Assign zones for this window
        std::vector<uint8_t> windowZones;
        assignZonesViterbi(window, view, windowZones, m_jumpTolerance, m_maxZoneCandidates);
        for (size_t i = 0; i < winLen; ++i)
            zones[start + i] = windowZones[i];
    }

    // Build a new PointView with processed points
    PointViewPtr newView = view->makeNew();
    PointRef srcPt = view->point(0);
    for (size_t shotIdx = 0; shotIdx < shotGroups.size(); ++shotIdx)
    {
        uint8_t zone = zones[shotIdx];
        auto& group = shotGroups[shotIdx];
        if (group.empty()) continue;
        std::vector<size_t> validEchoes;
        for (size_t idx : group)
        {
            srcPt.setPointId(idx);
            uint8_t nret = srcPt.getFieldAs<uint8_t>(Id::NumberOfReturns);
            if (nret > 0)
                validEchoes.push_back(idx);
        }
        if (validEchoes.empty())
            continue;
        // Reorder ReturnNumber and NumberOfReturns for valid echoes
        for (size_t i = 0; i < validEchoes.size(); ++i)
        {
            srcPt.setPointId(validEchoes[i]);
            srcPt.setField(Id::ReturnNumber, uint8_t(i + 1));
            srcPt.setField(Id::NumberOfReturns, uint8_t(validEchoes.size()));
        }
        // Apply zone correction and copy to newView
        for (size_t idx : validEchoes)
        {
            srcPt.setPointId(idx);
            PointId dstId = newView->size();
            newView->appendPoint(*view, idx);
            if (zone != 1)
            {
                double r = srcPt.getFieldAs<double>(Dimension::Id::EchoRange);
                double unambiguousRange = getUnambiguousRange(srcPt);
                double correction = (zone - 1) * unambiguousRange;
                double correctedRange = r + correction;
                double ox = srcPt.getFieldAs<double>(Dimension::Id::BeamOriginX);
                double oy = srcPt.getFieldAs<double>(Dimension::Id::BeamOriginY);
                double oz = srcPt.getFieldAs<double>(Dimension::Id::BeamOriginZ);
                double dx = srcPt.getFieldAs<double>(Dimension::Id::BeamDirectionX);
                double dy = srcPt.getFieldAs<double>(Dimension::Id::BeamDirectionY);
                double dz = srcPt.getFieldAs<double>(Dimension::Id::BeamDirectionZ);
                newView->setField(Dimension::Id::X, dstId, ox + dx * correctedRange);
                newView->setField(Dimension::Id::Y, dstId, oy + dy * correctedRange);
                newView->setField(Dimension::Id::Z, dstId, oz + dz * correctedRange);
            }
        }
    }
    PointViewSet viewSet;
    viewSet.insert(newView);
    return viewSet;
}

uint8_t MtaProcessingFilter::correctRangeGate(double measuredRange, double unambiguousRange)
{
        // Legacy code removed
    if (unambiguousRange <= 0) return 0; // Cannot determine without unambiguous range
    
    for (uint8_t zone = 1; zone <= m_maxZoneCandidates; ++zone)
    {
        double correction = (zone - 1) * unambiguousRange;
        double correctedRange = measuredRange + correction;
        
        if (correctedRange >= m_gateNear && correctedRange <= m_gateFar)
        {
            return zone;
        }
    }
    
    return 0; // No zone found within gate
}

double MtaProcessingFilter::getUnambiguousRange(PointRef& point) const
{
    // Priority:
    // 1. User-specified unambiguous_range (direct override)
    if (m_unambiguousRangeOverride > 0)
        return m_unambiguousRangeOverride;
    
    // 2. Read from dimension
    return point.getFieldAs<double>(Dimension::Id::UnambiguousRange);
}

void MtaProcessingFilter::assignZonesViterbi(const std::vector<ShotData>& window,
                                             PointViewPtr view,
                                             std::vector<uint8_t>& zonesOut,
                                             double jumpTolerance,
                                             uint8_t maxZones) const
{
    const size_t N = window.size();
    zonesOut.assign(N, 1);
    if (N == 0)
        return;

    // Helper: measured range and unambiguous range for a shot (from first echo)

    // DP tables: log-costs for numerical stability
    const uint8_t Z = std::max<uint8_t>(1, maxZones);
    std::vector<std::vector<double>> dp(N, std::vector<double>(Z, std::numeric_limits<double>::infinity()));
    std::vector<std::vector<uint8_t>> back(N, std::vector<uint8_t>(Z, 1));

    // Use view as an argument
    auto emissionCost = [&](size_t i, uint8_t z) -> double {
        // Favor zones where the corrected reflectance is consistent with the distance
        if (window[i].echoes.empty())
            return 1.0; // high cost for empty shot
        PointRef pt = view->point(0);
        size_t idx = window[i].echoes[0];
        pt.setPointId(idx);
        double measuredRange = pt.getFieldAs<double>(Dimension::Id::EchoRange);
        double unambiguousRange = getUnambiguousRange(pt);
        double correctedRange = measuredRange + (z - 1) * unambiguousRange;
        double reflectance = pt.getFieldAs<double>(Dimension::Id::Reflectance);
        // Expected attenuation law: reflectance ~ 1/(range^2)
        double expected = 1.0 / (correctedRange * correctedRange + 1e-6); // avoid div/0
        double relError = std::abs(reflectance - expected) / (expected + 1e-6);
        // Combine this cost with a slight preference for lower zones
        return relError + 0.01 * (z - 1);
    };



    auto transitionCost = [&](size_t prevIdx, uint8_t zp, size_t curIdx, uint8_t zc) -> double {
        // Optionally, use unambiguous range from first echo in each shot if needed
        // For now, just penalize zone jumps
        return 0.1 * std::abs(int(zc) - int(zp));
    };

    // Initialization
    for (uint8_t z = 1; z <= Z; ++z)
    {
        dp[0][z - 1] = emissionCost(0, z);
        back[0][z - 1] = z;
    }

    // Forward pass
    for (size_t i = 1; i < N; ++i)
    {
        for (uint8_t zc = 1; zc <= Z; ++zc)
        {
            double best = std::numeric_limits<double>::infinity();
            uint8_t bestZp = 1;
            for (uint8_t zp = 1; zp <= Z; ++zp)
            {
                double cost = dp[i - 1][zp - 1] + transitionCost(i-1, zp, i, zc) + emissionCost(i, zc);
                if (cost < best)
                {
                    best = cost;
                    bestZp = zp;
                }
            }
            dp[i][zc - 1] = best;
            back[i][zc - 1] = bestZp;
        }
    }

    // Backward reconstruction
    uint8_t z = 1;
    {
        double best = std::numeric_limits<double>::infinity();
        for (uint8_t zz = 1; zz <= Z; ++zz)
        {
            if (dp[N - 1][zz - 1] < best)
            {
                best = dp[N - 1][zz - 1];
                z = zz;
            }
        }
    }
    zonesOut[N - 1] = z;
    for (size_t i = N - 1; i > 0; --i)
    {
        z = back[i][z - 1];
        zonesOut[i - 1] = z;
    }
}

} // namespace pdal
