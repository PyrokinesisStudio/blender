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

/** \file gameengine/Ketsji/KX_SG_NodeRelationships.cpp
 *  \ingroup ketsji
 */


#include "KX_SG_NodeRelationships.h"

#include "BLI_utildefines.h"

/**
 * Implementation of classes defined in KX_SG_NodeRelationships.h
 */

/** 
 * first of all KX_NormalParentRelation
 */

	KX_NormalParentRelation *
KX_NormalParentRelation::
New(
) {
	return new KX_NormalParentRelation();
}

	bool
KX_NormalParentRelation::
UpdateChildCoordinates(
	SG_Node * child,
	const SG_Node * parent,
	bool& parentUpdated
) {
	BLI_assert(child != nullptr);

	if (!parentUpdated && !child->IsModified())
		return false;

	parentUpdated = true;

	if (parent==nullptr) { /* Simple case */
		child->SetWorldFromLocalTransform();
		child->ClearModified();
		return true; //false;
	}
	else {
		// the childs world locations which we will update.
		const mt::vec3 & p_world_scale = parent->GetWorldScaling();
		const mt::vec3 & p_world_pos = parent->GetWorldPosition();
		const mt::mat3 & p_world_rotation = parent->GetWorldOrientation();

		child->SetWorldScale(p_world_scale * child->GetLocalScale());
		child->SetWorldOrientation(p_world_rotation * child->GetLocalOrientation());
		child->SetWorldPosition(p_world_pos + p_world_scale * (p_world_rotation * child->GetLocalPosition()));
		child->ClearModified();
		return true;
	}
}

	SG_ParentRelation *
KX_NormalParentRelation::
NewCopy(
) {
	return new KX_NormalParentRelation();
}

KX_NormalParentRelation::
~KX_NormalParentRelation(
) {
	//nothing to do
}


KX_NormalParentRelation::
KX_NormalParentRelation(
) {
	// nothing to do
}

/** 
 * Next KX_VertexParentRelation
 */


	KX_VertexParentRelation *
KX_VertexParentRelation::
New(
) {
	return new KX_VertexParentRelation();
}
		
/** 
 * Method inherited from KX_ParentRelation
 */

	bool
KX_VertexParentRelation::
UpdateChildCoordinates(
	SG_Node * child,
	const SG_Node * parent,
	bool& parentUpdated
) {

	BLI_assert(child != nullptr);

	if (!parentUpdated && !child->IsModified())
		return false;

	child->SetWorldScale(child->GetLocalScale());
	
	if (parent)
		child->SetWorldPosition(child->GetLocalPosition()+parent->GetWorldPosition());
	else
		child->SetWorldPosition(child->GetLocalPosition());
	
	child->SetWorldOrientation(child->GetLocalOrientation());
	child->ClearModified();
	return true; //parent != nullptr;
}

/** 
 * Method inherited from KX_ParentRelation
 */

	SG_ParentRelation *
KX_VertexParentRelation::
NewCopy(
) {
	return new KX_VertexParentRelation();
};

KX_VertexParentRelation::
~KX_VertexParentRelation(
) {
	//nothing to do
}


KX_VertexParentRelation::
KX_VertexParentRelation(
) {
	//nothing to do
}


/**
 * Slow parent relationship
 */

	KX_SlowParentRelation *
KX_SlowParentRelation::
New(
	float relaxation
) {
	return new 	KX_SlowParentRelation(relaxation);
}

/** 
 * Method inherited from KX_ParentRelation
 */

	bool
KX_SlowParentRelation::
UpdateChildCoordinates(
	SG_Node * child,
	const SG_Node * parent,
	bool& parentUpdated
) {
	BLI_assert(child != nullptr);

	// the child will move even if the parent is not
	parentUpdated = true;

	const mt::vec3 & child_scale = child->GetLocalScale();
	const mt::vec3 & child_pos = child->GetLocalPosition();
	const mt::mat3 & child_rotation = child->GetLocalOrientation();

	// the childs world locations which we will update.
	
	mt::vec3 child_w_scale;
	mt::vec3 child_w_pos;
	mt::mat3 child_w_rotation;
		
	if (parent) {

		// This is a slow parent relation
		// first compute the normal child world coordinates.

		mt::vec3 child_n_scale;
		mt::vec3 child_n_pos;
		mt::mat3 child_n_rotation;

		const mt::vec3 & p_world_scale = parent->GetWorldScaling();
		const mt::vec3 & p_world_pos = parent->GetWorldPosition();
		const mt::mat3 & p_world_rotation = parent->GetWorldOrientation();

		child_n_scale = p_world_scale * child_scale;
		child_n_rotation = p_world_rotation * child_rotation;

		child_n_pos = p_world_pos + p_world_scale * 
			(p_world_rotation * child_pos);


		if (m_initialized) {

			// get the current world positions

			child_w_scale = child->GetWorldScaling();
			child_w_pos = child->GetWorldPosition();
			child_w_rotation = child->GetWorldOrientation();

			// now 'interpolate' the normal coordinates with the last 
			// world coordinates to get the new world coordinates.

			float weight = 1.0f / (m_relax + 1.0f);
			child_w_scale = (m_relax * child_w_scale + child_n_scale) * weight;
			child_w_pos = (m_relax * child_w_pos + child_n_pos) * weight;
			// for rotation we must go through quaternion
			const mt::quat child_w_quat = mt::quat::FromMatrix(child_w_rotation);
			const mt::quat child_n_quat = mt::quat::FromMatrix(child_n_rotation);
			child_w_rotation = mt::quat::Slerp(child_w_quat, child_n_quat, weight).ToMatrix();
			//FIXME: update physics controller.
		} else {
			child_w_scale = child_n_scale;
			child_w_pos = child_n_pos;
			child_w_rotation = child_n_rotation;
			m_initialized = true;
		}
			
	} else {

		child_w_scale = child_scale;
		child_w_pos = child_pos;
		child_w_rotation = child_rotation;
	}

	child->SetWorldScale(child_w_scale);
	child->SetWorldPosition(child_w_pos);
	child->SetWorldOrientation(child_w_rotation);
	child->ClearModified();
	// this node must always be updated, so reschedule it for next time
	child->ActivateRecheduleUpdateCallback();
	
	return true; //parent != nullptr;
}

/** 
 * Method inherited from KX_ParentRelation
 */

	SG_ParentRelation *
KX_SlowParentRelation::
NewCopy(
) {
	return new 	KX_SlowParentRelation(m_relax);
}

KX_SlowParentRelation::
KX_SlowParentRelation(
	float relaxation
):
	m_relax(relaxation),
	m_initialized(false)
{
	//nothing to do
}

KX_SlowParentRelation::
~KX_SlowParentRelation(
) {
	//nothing to do
}




