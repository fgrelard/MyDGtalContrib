#ifndef DIGITALPLANESET_H
#define DIGITALPLANESET_H

#include "DigitalPlane.h"

template <typename TSpace>
class DigitalPlaneSet {
public:
        typedef DigitalPlane<TSpace> Plane;
        typedef typename Plane::DigitalSet DigitalSet;
public:
        DigitalPlaneSet(const Plane& digitalPlane,
                        const DigitalSet& aDigitalSet) : myDigitalPlane(digitalPlane),
                                                       myDigitalSet(aDigitalSet) {}

        DigitalPlaneSet(const DigitalPlaneSet& other) : myDigitalPlane(other.myDigitalPlane), myDigitalSet(other.myDigitalSet) {}

        DigitalSet pointSet() const { return myDigitalSet; }
        Plane digitalPlane() const { return myDigitalPlane; }
private:
        Plane myDigitalPlane;
        DigitalSet myDigitalSet;
};

#endif
