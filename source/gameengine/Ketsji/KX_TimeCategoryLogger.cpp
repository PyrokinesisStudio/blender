/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file gameengine/Ketsji/KX_TimeCategoryLogger.cpp
 *  \ingroup ketsji
 */

#include "KX_TimeCategoryLogger.h"
#include "RAS_DebugDraw.h"

#include "boost/format.hpp"

const std::string profileLabels[KX_TimeLogger::NUM_CATEGORY] = {
	"Physics", // PHYSICS
	"Logic", // LOGIC
	"Animation", // ANIMATIONS
	"Network", // NETWORK
	"Scenegraph", // SCENEGRAPH
	"Rasterizer", // RASTERIZER
	"Services", // SERVICES
	"Overhead", // OVERHEAD
	"Outside", // OUTSIDE
	"GPU Latency" // LATENCY
};

KX_TimeCategoryLogger::KX_TimeCategoryLogger()
	:m_lastCategory(KX_TimeLogger::NONE),
	m_lastTotalAverage(0.0)
{
	for (double& avg : m_lastAverages) {
		avg = 0.0;
	}
}

KX_TimeCategoryLogger::~KX_TimeCategoryLogger()
{
}

void KX_TimeCategoryLogger::StartLog(KX_TimeLogger::Category tc, double now)
{
	if (m_lastCategory != KX_TimeLogger::NONE) {
		m_loggers[m_lastCategory].EndLog(now);
	}
	m_loggers[tc].StartLog(now);
	m_lastCategory = tc;
}

void KX_TimeCategoryLogger::EndLog(KX_TimeLogger::Category tc, double now)
{
	m_loggers[tc].EndLog(now);
}

void KX_TimeCategoryLogger::EndLog(double now)
{
	m_loggers[m_lastCategory].EndLog(now);
	m_lastCategory = KX_TimeLogger::NONE;
}

void KX_TimeCategoryLogger::NextMeasurement(double now)
{
	m_lastTotalAverage = 0.0;
	for (unsigned short tc = 0; tc < KX_TimeLogger::NUM_CATEGORY; ++tc) {
		KX_TimeLogger& logger = m_loggers[tc];
		logger.NextMeasurement(now);

		const double time = logger.GetAverage();
		m_lastAverages[tc] = time;
		m_lastTotalAverage += time;
	}
}

double KX_TimeCategoryLogger::GetAverage(KX_TimeLogger::Category tc) const
{
	return m_lastAverages[tc];
}

double KX_TimeCategoryLogger::GetAverage() const
{
	return m_lastTotalAverage;
}

double KX_TimeCategoryLogger::GetAverageFrameRate() const
{
	if (m_lastTotalAverage < 1e-6) {
		// Equivalent to 1.0 / 1e-6.
		return 1e6f;
	}
	return 1.0 / m_lastTotalAverage;
}

std::map<std::string, double> KX_TimeCategoryLogger::GetProfileDict()
{
	std::map<std::string, double> dict;

	for (unsigned short tc = 0; tc < KX_TimeLogger::NUM_CATEGORY; ++tc) {
		dict[profileLabels[tc]] = m_lastAverages[tc];
	}

	return dict;
}

static const MT_Vector4 white(1.0f, 1.0f, 1.0f, 1.0f);

void KX_TimeCategoryLogger::RenderFrameRate(RAS_DebugDraw& debugDraw, int xindent, int ysize,
											 int& xcoord, int& ycoord, int profileIndent)
{
	debugDraw.RenderText2D("Frametime :", MT_Vector2(xcoord + xindent, ycoord), white);

	const std::string debugtxt = (boost::format("%5.2fms (%.1ffps)") %  (m_lastTotalAverage * 1000.0f) % GetAverageFrameRate()).str();
	debugDraw.RenderText2D(debugtxt, MT_Vector2(xcoord + xindent + profileIndent, ycoord), white);
	// Increase the indent by default increase
	ycoord += ysize;
}

void KX_TimeCategoryLogger::RenderCategories(RAS_DebugDraw& debugDraw, int xindent, int ysize,
											 int& xcoord, int& ycoord, int profileIndent)
{
	double tottime = m_lastTotalAverage;
	if (tottime < 1e-6) {
		tottime = 1e-6;
	}

	for (unsigned short tc = 0; tc < KX_TimeLogger::NUM_CATEGORY; ++tc) {
		debugDraw.RenderText2D(profileLabels[tc] + ":", MT_Vector2(xcoord + xindent, ycoord), white);

		const double time = m_lastAverages[tc];

		const std::string debugtxt = (boost::format("%5.2fms | %d%%") % (time*1000.f) % (int)(time/tottime * 100.f)).str();
		debugDraw.RenderText2D(debugtxt, MT_Vector2(xcoord + xindent + profileIndent, ycoord), white);

		const MT_Vector2 boxSize(50 * (time / tottime), 10);
		debugDraw.RenderBox2D(MT_Vector2(xcoord + (int)(2.2 * profileIndent), ycoord), boxSize, white);
		ycoord += ysize;
	}
}
