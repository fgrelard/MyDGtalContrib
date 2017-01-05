#ifndef CURVE_H
#define CURVE_H

template <typename Container>
class Curve {
public:
        typedef typename Container::value_type Point;
public:
        Curve() = delete;
        Curve(const Container& aCurve) : myPoints(aCurve) {}
        Curve(const Curve<Container>& other) : myPoints(other.myPoints) {}
        Container pointSet() const { return myPoints; }
private:
        Container myPoints;
};

#endif
