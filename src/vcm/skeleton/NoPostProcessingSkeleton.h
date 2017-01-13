#ifndef NO_POST_PROCESSING_SKELETON_H
#define NO_POST_PROCESSING_SKELETON_H

template <typename Container>
class NoPostProcessingSkeleton {
public:
        NoPostProcessingSkeleton() = delete;
        NoPostProcessingSkeleton(const Container& skeletonPoints) {
                mySkeleton = new Container( skeletonPoints);
        }
        NoPostProcessingSkeleton(const NoPostProcessingSkeleton& other) {
                mySkeleton = new Container (*other.mySkeleton);
        }
        ~NoPostProcessingSkeleton() {
                if (mySkeleton != 0) {
                        delete mySkeleton;
                        mySkeleton = 0;
                }
        }

public:
        Container postProcess() { return *mySkeleton; }

private:
        Container* mySkeleton;
};

#endif
