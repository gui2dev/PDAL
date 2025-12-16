/*
 * Unit tests for pdal::MtaProcessingFilter
 *
 * Copyright (c) 2025
 */

#include <pdal/pdal_test_main.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/Options.hpp>
#include <pdal/Filter.hpp>
#include <pdal/filters/MtaProcessingFilter.hpp>
#include <gtest/gtest.h>

using namespace pdal;

namespace {

// Helper to create a simple PointView with required dimensions
PointViewPtr createTestView(PointTable& table, const std::vector<double>& ranges, double ua, double reflectance = 1.0) {
    table.layout()->registerDim(Dimension::Id::X);
    table.layout()->registerDim(Dimension::Id::Y);
    table.layout()->registerDim(Dimension::Id::Z);
    table.layout()->registerDim(Dimension::Id::EchoRange);
    table.layout()->registerDim(Dimension::Id::UnambiguousRange);
    table.layout()->registerDim(Dimension::Id::Reflectance);
    table.layout()->registerDim(Dimension::Id::BeamOriginX);
    table.layout()->registerDim(Dimension::Id::BeamOriginY);
    table.layout()->registerDim(Dimension::Id::BeamOriginZ);
    table.layout()->registerDim(Dimension::Id::BeamDirectionX);
    table.layout()->registerDim(Dimension::Id::BeamDirectionY);
    table.layout()->registerDim(Dimension::Id::BeamDirectionZ);
    table.layout()->registerDim(Dimension::Id::GpsTime);
    table.layout()->registerDim(Dimension::Id::ReturnNumber);
    table.layout()->registerDim(Dimension::Id::NumberOfReturns);

    PointViewPtr view(new PointView(table));
    for (size_t i = 0; i < ranges.size(); ++i) {
        PointId id = view->size();
        view->appendPoint(*view, id); // dummy, will be overwritten
        view->setField(Dimension::Id::EchoRange, id, ranges[i]);
        view->setField(Dimension::Id::UnambiguousRange, id, ua);
        view->setField(Dimension::Id::Reflectance, id, reflectance);
        view->setField(Dimension::Id::BeamOriginX, id, 0.0);
        view->setField(Dimension::Id::BeamOriginY, id, 0.0);
        view->setField(Dimension::Id::BeamOriginZ, id, 0.0);
        view->setField(Dimension::Id::BeamDirectionX, id, 1.0);
        view->setField(Dimension::Id::BeamDirectionY, id, 0.0);
        view->setField(Dimension::Id::BeamDirectionZ, id, 0.0);
        view->setField(Dimension::Id::GpsTime, id, static_cast<double>(i));
        view->setField(Dimension::Id::ReturnNumber, id, 1);
        view->setField(Dimension::Id::NumberOfReturns, id, 1);
    }
    return view;
}

} // namespace

TEST(MtaProcessingFilterTest, FixedMode)
{
    PointTable table;
    auto view = createTestView(table, {10.0, 20.0, 30.0}, 15.0);
    Options opts;
    opts.add("mode", "fixed");
    opts.add("fixed_zone", 2);
    MtaProcessingFilter filter;
    filter.setOptions(opts);
    filter.prepare(table);
    PointViewSet out = filter.run(view);
    ASSERT_EQ(out.size(), 1U);
    auto outView = *out.begin();
    ASSERT_EQ(outView->size(), 3U);
    // All points should be shifted by one UA (zone 2)
    for (PointId i = 0; i < 3; ++i) {
        double x = outView->getFieldAs<double>(Dimension::Id::X, i);
        EXPECT_DOUBLE_EQ(x, 10.0 + i * 10.0 + 15.0); // original + UA
    }
}

TEST(MtaProcessingFilterTest, AutoModeViterbi)
{
    PointTable table;
    // Simulate ambiguous ranges: first two in zone 1, last in zone 2
    auto view = createTestView(table, {10.0, 12.0, 25.0}, 15.0, 1.0);
    Options opts;
    opts.add("mode", "auto");
    opts.add("max_zones", 2);
    MtaProcessingFilter filter;
    filter.setOptions(opts);
    filter.prepare(table);
    PointViewSet out = filter.run(view);
    ASSERT_EQ(out.size(), 1U);
    auto outView = *out.begin();
    ASSERT_EQ(outView->size(), 3U);
    // The last point should be assigned to zone 2 (shifted by UA)
    double x0 = outView->getFieldAs<double>(Dimension::Id::X, 0);
    double x1 = outView->getFieldAs<double>(Dimension::Id::X, 1);
    double x2 = outView->getFieldAs<double>(Dimension::Id::X, 2);
    EXPECT_NEAR(x0, 10.0, 1e-6);
    EXPECT_NEAR(x1, 12.0, 1e-6);
    EXPECT_NEAR(x2, 25.0 + 15.0, 1e-6);
}

TEST(MtaProcessingFilterTest, RangeGateMode)
{
    PointTable table;
    // Points: one in gate, one below, one above
    auto view = createTestView(table, {5.0, 15.0, 25.0}, 10.0);
    Options opts;
    opts.add("mode", "rangegate");
    opts.add("gate_near", 10.0);
    opts.add("gate_far", 20.0);
    opts.add("max_zones", 3);
    MtaProcessingFilter filter;
    filter.setOptions(opts);
    filter.prepare(table);
    PointViewSet out = filter.run(view);
    ASSERT_EQ(out.size(), 1U);
    auto outView = *out.begin();
    ASSERT_EQ(outView->size(), 3U);
    // Only the second point should be in the gate (zone 1, 15.0)
    // The first and third should be outside, so zone 0 (no correction)
    double x0 = outView->getFieldAs<double>(Dimension::Id::X, 0);
    double x1 = outView->getFieldAs<double>(Dimension::Id::X, 1);
    double x2 = outView->getFieldAs<double>(Dimension::Id::X, 2);
    // x1 should be 15.0 (in gate), x0 and x2 unchanged
    EXPECT_DOUBLE_EQ(x0, 5.0);
    EXPECT_DOUBLE_EQ(x1, 15.0);
    EXPECT_DOUBLE_EQ(x2, 25.0);
}
