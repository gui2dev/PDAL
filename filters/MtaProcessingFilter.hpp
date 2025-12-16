/******************************************************************************
 * Copyright (c) 2025, Guilhem Villemin
 * All rights reserved.
 ****************************************************************************/

#pragma once

#include <pdal/Filter.hpp>
#include <pdal/Streamable.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/util/ProgramArgs.hpp>
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace pdal
{

struct ShotData;


class PDAL_EXPORT MtaProcessingFilter : public Filter, public Streamable
{
public:
    MtaProcessingFilter();
    ~MtaProcessingFilter() override;

    std::string getName() const override;

private:
    void addArgs(ProgramArgs& args) override;
    void initialize() override;
    void prepared(PointTableRef table) override;
    bool processOne(PointRef& point) override;
    PointViewSet run(PointViewPtr view) override;

    uint8_t correctRangeGate(double measuredRange, double unambiguousRange);
    double getUnambiguousRange(PointRef& point) const;
    void assignZonesViterbi(const std::vector<ShotData>& window,
                            PointViewPtr view,
                            std::vector<uint8_t>& zonesOut,
                            double jumpTolerance,
                            uint8_t maxZones) const;



    std::string m_mode;
    uint8_t m_fixedZone;
    float m_gateNear;
    float m_gateFar;
    uint8_t m_maxZoneCandidates;
    double m_unambiguousRangeOverride;
    double m_maxOriginDistance;
    double m_maxAngleDifference;
    double m_jumpTolerance;
};

} // namespace pdal
