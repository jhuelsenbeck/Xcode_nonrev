#ifndef TransitionProbabilities_H
#define TransitionProbabilities_H

class Model;


class TransitionProbabilities {

    public:
                                    TransitionProbabilities(void) = delete;
                                    TransitionProbabilities(int s, int k, Model* m);
                                   ~TransitionProbabilities(void);
        TransitionProbabilities&    operator=(TransitionProbabilities& c);
        void                        print(void);
        void                        tiProbs(double v);
        double***                   getTiProbs(void) { return p; }
    
    protected:
        Model*                      myModel;
        double***                   p;
        int                         mySpace;
        int                         numStates;
        int                         numGammaCats;
        int                         tiSizeInBytes;
};

#endif
