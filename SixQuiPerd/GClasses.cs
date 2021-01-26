using System;
using System.Windows.Forms;
using System.Drawing;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace SixQuiPerd
{
    public static class Extentions
    {
        public static void Each<T>(this IEnumerable<T> items, Action<T> action)
        {
            foreach (var item in items)
            {
                action(item);
            }
        }

        public static IEnumerable<T> Select<T>(int n, Func<T> function)
        {
            T[] _tab = new T[n];
            for (int i = 0; i < n; ++i)
                _tab[i] = function.Invoke();
            return _tab;
        }
    }

    struct CInterval
    {
        public int a, b;
    }
    interface IGenJr
    {
        void setJeu(CPlateau plat, CPioche ppioch, int pjrNbr);
        void reInitialiser();
        byte programmer();
        int choisirLigne(byte cr);
        void recuperer(byte c);
        void retirer(byte c);
        void jrAction(int jr, byte cr, int ln);
        void piocheMelange();
        void piocheMelange(byte cr);
        int nbCrdMain();
        void scoreDelta(int nmJr, int[] delta);
    }

    class CGenJ : IGenJr
    {
        protected Random rnd;
        protected List<byte> cds;
        protected CPlateau plateau;
        protected CPioche pioch;
        protected int jrNbr;

        public CGenJ()
        {
            cds = new List<byte>();
            plateau = null;
            jrNbr = 0;
            rnd = new Random();
        }

        virtual public void reInitialiser()
        {
            cds.Clear();
            rnd = new Random(rnd.Next()^(new Random().Next()));
        }

        virtual public void recuperer(byte c){cds.Add(c);}

        virtual public int choisirLigne(byte cr)
        {
            //Random rnd = new Random(DateTime.Now.Millisecond);
            int res = (int)(Math.Floor(4 * rnd.NextDouble()));
            if (res >= 4) res = 4 - 1;
            return res;
        }

        virtual public void piocheMelange() { }
        virtual public void piocheMelange(byte cr) { }

        virtual public void jrAction(int jr, byte cr, int ln){}

        virtual public byte programmer()
        {
            if (cds.Count > 0)
            {
                //Random rnd = new Random(DateTime.Now.Millisecond);
                int num = (int)(Math.Floor(cds.Count*rnd.NextDouble()));
                if (num >= cds.Count) num = cds.Count - 1;
                foreach (byte cr in cds)
                {
                    if (1 <= cr && cr <= 104)
                    {
                        if (num == 0) return cr;
                        else --num;
                    }
                }
            }
            return 0;
        }

        virtual public void retirer(byte c){cds.Remove(c);}
        virtual public void setJeu(CPlateau plat, CPioche ppioch, int pjrNbr) { plateau = plat; pioch = ppioch; jrNbr = pjrNbr; }
        virtual public int nbCrdMain() { return cds.Count; }

        virtual public void scoreDelta(int nmJr, int[] delta){}
    }
    class CMalusJ : CGenJ
    {
        override public int choisirLigne(byte cr)
        {
            int res = 0;
            int min = plateau.malus(0);
            for (int l = 1; l < 4; ++l)
            {
                int m = plateau.malus(l);
                if (m< min)
                {
                    min = m;
                    res = l;
                }
            }
            return res;
        }
    }
    class CProgMalusJ : CMalusJ
    {
        double pondPrend;

        public CProgMalusJ(double ppondPrend = 1.0f)
        {
            pondPrend = ppondPrend;
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                double min = double.MaxValue;
                byte res = 0;

                int minremp = int.MaxValue;
                int[] malus = new int[4];
                for (int l = 0; l < 4; ++l)
                {
                    malus[l] = plateau.malus(l);
                    minremp = Math.Min(minremp, malus[l]);
                }

                foreach (byte cr in cds)
                {
                    if (1 <= cr && cr <= 104)
                    {
                        int l = plateau.ouPlacer(cr);
                        double malu = (l >= 0 ? malus[l] : minremp * pondPrend);
                        if (malu < min)
                        {
                            res = cr;
                            min = malu;
                        }
                    }
                }
                return res;
            }
            return 0;
        }
    }
    class CProgJ : CMalusJ
    {
        double pondPrend;
        double pondRisq6P;

        public CProgJ(double ppondPrend = 1.0f, double ppondRisq6P = 1.0f)
        {
            pondPrend = ppondPrend;
            pondRisq6P = ppondRisq6P;
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                double min = double.MaxValue;
                byte res = 0;

                int minremp = int.MaxValue;
                int[] malus = new int[4];
                int[] derCrd = new int[4];
                int[] nbCrd = new int[4];
                for (int l = 0; l < 4; ++l)
                {
                    malus[l] = plateau.malus(l);
                    minremp = Math.Min(minremp, malus[l]);
                    derCrd[l] = plateau.derniereCrd(l);
                    nbCrd[l] = plateau.nbCartes(l);
                }

                foreach (byte cr in cds)
                {
                    if (1 <= cr && cr <= 104)
                    {
                        int l = plateau.ouPlacer(cr);
                        double malu;
                        if (l >= 0)
                        {
                            double maluSaut; //Estimation du risque si la ligne saute
                            if (malus[l] == minremp)//si il est possible que la ligne saute ?
                            {
                                int maxLb = 0;
                                int lb = -1;
                                for (int i = 0; i < 4; ++i) if (i != l)
                                    {
                                        if (maxLb < derCrd[i] && derCrd[i] < cr)
                                        {
                                            maxLb = derCrd[i];
                                            lb = i;
                                        }
                                    }
                                if (lb >= 0)
                                {
                                    if (nbCrd[lb] + (Math.Min(cr - derCrd[lb], plateau.obtenirNbJ() - 1)) < 6)//On ne rammassera pas !
                                        maluSaut = 0.0f;
                                    else maluSaut = malus[l];
                                }
                                else maluSaut = 0.0f;
                            }
                            else maluSaut = 0.0f;

                            if (nbCrd[l] + (Math.Min(cr - derCrd[l], plateau.obtenirNbJ())) < 6)//On ne rammassera pas !
                                malu = 0.0f;
                            else malu = malus[l];

                            malu += maluSaut * pondRisq6P;
                        }
                        else malu = minremp * pondPrend;
                        if (malu < min)
                        {
                            res = cr;
                            min = malu;
                        }
                    }
                }
                return res;
            }
            return 0;
        }
    }
    class CProgCptJ : CMalusJ
    {
        protected double pondPrend;
        protected double pondRisq6P;

        protected byte[] crdsMem;

        public CProgCptJ(double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
        {
            pondPrend = ppondPrend;
            pondRisq6P = ppondRisq6P;
            crdsMem = new byte[104];
            for (int i = 0; i < 104; ++i) crdsMem[i] = 255;
        }

        override public void recuperer(byte c)
        {
            base.recuperer(c);
            if (1 <= c && c <= 104) crdsMem[c-1] = 0;
        }

        override public void reInitialiser()
        {
            base.reInitialiser();
            for (int i = 0; i < 104; ++i) crdsMem[i] = 255;
        }

        override public void jrAction(int jr, byte cr, int ln)
        {
            if(1<=cr && cr<=104) crdsMem[cr-1] = 0;
        }

        override public void piocheMelange() { for (int i = 0; i < 104; ++i) crdsMem[i] = 255; }
        override public void piocheMelange(byte cr) { if(1<=cr && cr<=104) crdsMem[cr-1] = 255; }

        protected int nbCrdRstItr(byte ca, byte cb)
        {
            int res = 0;
            for (int i = ca /*+ 1 - 1*/; i < ((int)cb)-1; ++i) if (crdsMem[i] == 255) ++res;
            return res;
        }

        /*protected double calculProbaCrdAutreJr(byte ca, byte cb)
        {
            int nbCardJr = cds.Count * (plateau.obtenirNbJ() - 1); //Nombre total de cartes parmis les autres joueurs
            int nbCrdTot = nbCardJr + pioch.NbCartes(); //Nombre total de cartes inconnues

            if (nbCrdTot == 0) return 0;
            return (nbCrdRstItr(ca,cb)*nbCardJr) / (double)nbCrdTot;
        }*/

        protected double calculProbaCrdJr(int nbCrdRst, int nbJ)
        {
            int nbCardJr = cds.Count * nbJ; //Nombre total de cartes parmis les autres joueurs
            int nbCrdTot = nbCardJr + pioch.NbCartes(); //Nombre total de cartes inconnues

            if (nbCrdTot == 0) return 0.0f;
            return (nbCrdRst * nbCardJr) / (double)nbCrdTot;
        }

        protected void preCalculRsq(ref double minremp, int[] malus, byte[] derCrd, ref byte derCrdMin, int[] nbCrd, CPlateau plat)
        {
            minremp = double.MaxValue;
            derCrdMin = byte.MaxValue;
            for (int l = 0; l < 4; ++l)
            {
                malus[l] = plat.malus(l);
                derCrd[l] = plat.derniereCrd(l);
                minremp = Math.Min(minremp, (double)malus[l]);
                derCrdMin = Math.Min(derCrdMin, derCrd[l]);
                nbCrd[l] = plat.nbCartes(l);
            }
        }

        virtual protected double estimRisq(byte cr, double minremp, int[] malus, byte[] derCrd, byte derCrdMin, int[] nbCrd, CPlateau plat)
        {
            int l = plat.ouPlacer(cr);
            double malu;
            if (l >= 0)
            {
                double maluSaut; //Estimation du risque si la ligne saute
                                 // nbCrdRstItr(ca, cb) = Nombre de cartes restantes entre ca exclu et cb exclu
                if (malus[l] == minremp && nbCrdRstItr(1, derCrdMin) > 0)//si il est possible que la ligne saute ?
                {
                    int maxLb = 0;
                    int lb = -1;
                    for (int i = 0; i < 4; ++i) if (i != l)
                        {
                            if (maxLb < derCrd[i] && derCrd[i] < cr)
                            {
                                maxLb = derCrd[i];
                                lb = i;
                            }
                        }
                    if (lb >= 0)
                    {
                        // nbCrdRstItr(ca, cb) = Nombre de cartes restantes entre ca exclu et cb exclu
                        if (nbCrd[lb] + (Math.Min(1 + nbCrdRstItr(derCrd[lb], cr), plat.obtenirNbJ() - 1)) < 6)//On ne rammassera pas !
                            maluSaut = 0.0;
                        else maluSaut = malus[l];
                    }
                    else maluSaut = 0.0;
                }
                else maluSaut = 0.0f;
                // nbCrdRstItr(ca, cb) = Nombre de cartes restantes entre ca exclu et cb exclu
                if (nbCrd[l] + (Math.Min(1 + nbCrdRstItr(derCrd[l], cr), plat.obtenirNbJ())) < 6)//On ne rammassera pas !
                    malu = 0.0;
                else malu = malus[l];

                malu += maluSaut * pondRisq6P;
            }
            else malu = minremp * pondPrend;
            return malu;
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                double min = double.MaxValue;
                byte res = 0;

                double minremp = double.MaxValue;
                int[] malus = new int[4];
                byte[] derCrd = new byte[4];
                byte derCrdMin = byte.MaxValue;
                int[] nbCrd = new int[4];

                preCalculRsq(ref minremp,malus,derCrd,ref derCrdMin,nbCrd,plateau);

                foreach (byte cr in cds)
                {
                    if (1<= cr && cr<=104)
                    {
                        double malu = estimRisq(cr, minremp, malus, derCrd, derCrdMin, nbCrd, plateau);
                        if (malu < min)
                        {
                            res = cr;
                            min = malu;
                        }
                    }
                }
                return res;
            }
            return 0;
        }

    }

    class CProgCptJRun : CProgCptJ
    {
        public CProgCptJRun(double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
            :base(ppondPrend, ppondRisq6P)
        {
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                double min = double.MaxValue;
                byte crSup = 0;
                byte crInf = 0;

                double minremp = double.MaxValue;
                int[] malus = new int[4];
                byte[] derCrd = new byte[4];
                byte derCrdMin = byte.MaxValue;
                int[] nbCrd = new int[4];

                preCalculRsq(ref minremp, malus, derCrd, ref derCrdMin, nbCrd, plateau);

                foreach (byte cr in cds)
                {
                    if (1 <= cr && cr <= 104)
                    {
                        int l = plateau.ouPlacer(cr);
                        if(l>=0)
                        {
                            double malu = estimRisq(cr, minremp, malus, derCrd, derCrdMin, nbCrd, plateau);
                            if (malu < min)
                            {
                                crSup = cr;
                                min = malu;
                            }
                            else if(malu == 0.0)
                            {
                                crSup = Math.Min(crSup, cr);
                            }
                        }
                        else if(crInf < cr)
                        {
                            crInf = cr;
                        }
                    }
                }
                byte res = 0;

                if (crInf > 0)
                {
                    double malusInf = minremp * pondPrend;
                    if(malusInf < min)
                    {
                        res = crInf;
                    }
                    else
                    {
                        res = crSup;
                    }
                }
                else res = crSup;

                return res;
            }
            return 0;
        }
    }
    class CProgCptJRunSup : CProgCptJ
    {
        public CProgCptJRunSup(double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
            : base(ppondPrend, ppondRisq6P)
        {
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                int nbItrvfCrds = 0;
                double min = double.MaxValue;
                byte crSup = 0;
                byte crInf = 0;

                double minremp = double.MaxValue;
                int[] malus = new int[4];
                byte[] derCrd = new byte[4];
                byte derCrdMin = byte.MaxValue;
                int[] nbCrd = new int[4];

                preCalculRsq(ref minremp, malus, derCrd, ref derCrdMin, nbCrd, plateau);

                //byte[] choixCrLn = new byte[4];

                foreach (byte cr in cds)
                {
                    if (1 <= cr && cr <= 104)
                    {
                        int l = plateau.ouPlacer(cr);
                        if (l >= 0)
                        {
                            int tNbItrvfCrds = nbCrdRstItr(derCrd[l], cr);
                            double malu = estimRisq(cr, minremp, malus, derCrd, derCrdMin, nbCrd, plateau);
                            if (malu < min)
                            {
                                crSup = cr;
                                min = malu;
                                //choixCrLn[l] = cr;
                                nbItrvfCrds = tNbItrvfCrds;
                            }
                            else if (malu == min)
                            {
                                if(malu == 0.0)
                                {
                                    crSup = Math.Min(crSup, cr);
                                    nbItrvfCrds = 0;
                                }
                                else if (tNbItrvfCrds > nbItrvfCrds)
                                {
                                    crSup = cr;
                                    min = malu;
                                    //choixCrLn[l] = cr;
                                    nbItrvfCrds = tNbItrvfCrds;
                                }
                            }
                        }
                        else if (crInf < cr)
                        {
                            crInf = cr;
                        }
                    }
                }
                byte res = 0;

                if (crInf > 0)
                {
                    double malusInf = minremp * pondPrend;
                    if (malusInf < min)
                    {
                        res = crInf;
                    }
                    else
                    {
                        res = crSup;
                    }
                }
                else res = crSup;

                return res;
            }
            return 0;
        }
    }

    class CProgRecurJ : CProgCptJ
    {
        private double coefEstp;
        private int lig;
        private int depth;

        public CProgRecurJ(double ppondPrend = 1.0f, double ppondRisq6P = 1.0, int pdepth = 5, double pcoefEstp = 1.0)
            : base(ppondPrend, ppondRisq6P) { depth = pdepth; coefEstp = pcoefEstp; }


        protected byte[] duplic(byte[] orig)
        {
            byte[] res = new byte[orig.Length];
            for (int i = 0; i < orig.Length; ++i) res[i] = orig[i];
            return res;
        }

        protected List<byte> duplic(List<byte> orig)
        {
            List<byte> res = new List<byte>();
            foreach (byte b in orig) res.Add(b);
            return res;
        }

        private double estimRsqRecur(CPlateau plat, ref byte minCr, ref int lig, int dpth)
        {
            if (cds.Count > 0 && dpth>0)
            {
                minCr = 0;
                double min = double.MaxValue;
                int delta;

                double minremp = double.MaxValue;
                int[] malus = new int[4];
                byte[] derCrd = new byte[4];
                byte derCrdMin = byte.MaxValue;
                int[] nbCrd = new int[4];

                preCalculRsq(ref minremp, malus, derCrd, ref derCrdMin, nbCrd, plateau);

                List<byte> svCds = cds;
                foreach (byte cr in svCds)
                {
                    double malu = estimRisq(cr, minremp, malus, derCrd, derCrdMin, nbCrd, plateau);
                    //byte[] svCrdsMem = crdsMem;

                    cds = duplic(svCds);
                    cds.Remove(cr);
                    //crdsMem = duplic(svCrdsMem);

                    int l = plat.ouPlacer(cr);
                    if (l >= 0)
                    {
                        lig = l;
                        CPlateau fant = plat.obtenirFantome();
                        fant.placer(jrNbr, cr, -1, out delta);
                        byte nCrMin = 0;
                        int nLg = -1;
                        malu += coefEstp * estimRsqRecur(fant, ref nCrMin, ref nLg, dpth-1);
                        if (malu < min)
                        {
                            min = malu;
                            minCr = cr;
                        }
                    }
                    else
                    {
                        //double maluOrig = malu;
                        for(l = 0;l<4;++l)
                        {
                            malu = malus[l] * pondPrend;
                            CPlateau fant = plat.obtenirFantome();
                            fant.placer(jrNbr, cr, l, out delta);
                            byte nCrMin = 0;
                            int nLg = -1;
                            malu += coefEstp * estimRsqRecur(fant, ref nCrMin, ref nLg, dpth-1);
                            if (malu < min)
                            {
                                min = malu;
                                minCr = cr;
                                lig = l;
                            }
                        }
                    }
                    //crdsMem = svCrdsMem;
                }
                cds = svCds;

                return min;
            }
            else return 0.0;
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                byte cr = 0;
                lig = -1;
                double min = estimRsqRecur(plateau, ref cr, ref lig, depth);
                return cr;
            }
            else return 0;
        }

        override public int choisirLigne(byte cr) { return lig; }

    }
    class CProgProbJ : CProgCptJ
    {
        private double seuilSaut;
        private double puissSaut;

        public CProgProbJ(double ppondPrend = 1.0f, double ppondRisq6P = 1.0, double pseuilSaut = 0.9, double ppuissSaut = 0.05)
            : base(ppondPrend, ppondRisq6P){ seuilSaut = pseuilSaut; puissSaut = ppuissSaut; }

        private double probItrvl(double p, int nbCrds, double seuilInf = 1.0, double seuilSup = 1.0, double seuilSupMax = 10.0, double puissInf = 0.01f, double puissSup = 0.01f)
        {
            //seuilSupMax *= plateau.obtenirNbJ();
            p = p * (plateau.obtenirNbJ() - 1);
            //p /= (5 - nbCrds);
            nbCrds = 5 - nbCrds;
            seuilInf *= nbCrds;
            seuilSup += nbCrds;
            seuilSupMax += nbCrds;

            if (seuilInf <= p && p <= seuilSup) p = 1.0f;
            else if (seuilInf < p) p = Math.Pow( p / seuilInf, puissInf);
            else if (p >= seuilSupMax) p = 0.0f;
            else /*if (seuilSup < p && p < seuilSupMax)*/ p = Math.Pow((seuilSupMax-p) / (seuilSupMax - seuilSup), puissSup);
            //return 0.5+0.5*(1.0-p);

            return p;

            /*if (5.5f <= p && p <= 6.5f) p = 1.0f;
            else if (p > (6.0f + 5.5f)) p = 0.0f;
            else if (5.5f < p) p = p / 5.5f;
            else if (p > 6.5f) p = (6.5f-p) / 5.5f;
            return p;*/
        }

        private double probSaut(double p, double seuil = 0.9, double puiss = 0.05)
        {
            //return 1.0f;
            p = p * (plateau.obtenirNbJ() - 1);
            if (p >= seuil) p = 1.0f;
            else p = (double)Math.Pow(p/seuil, puiss);
            return p;
        }

        private double probFct(double p)
        {
            return p;
        }

        override public byte programmer()
        {
            if (cds.Count > 0)
            {
                double min = double.MaxValue;
                byte res = 0;

                double minremp = double.MaxValue;
                int[] malus = new int[4];
                byte[] derCrd = new byte[4];
                byte derCrdMin = byte.MaxValue;
                int[] nbCrd = new int[4];
                for (int l = 0; l < 4; ++l)
                {
                    malus[l] = plateau.malus(l);
                    derCrd[l] = plateau.derniereCrd(l);
                    minremp = Math.Min(minremp, (double)malus[l]);
                    derCrdMin = Math.Min(derCrdMin, derCrd[l]);
                    nbCrd[l] = plateau.nbCartes(l);
                }

                int nbCrdRst;
                foreach (byte cr in cds)
                {
                    int l = plateau.ouPlacer(cr);
                    double malu;
                    if (l >= 0)
                    {
                        double maluSaut; //Estimation du risque si la ligne saute
                        // nbCrdRstItr(ca, cb) = Nombre de cartes restantes entre ca exclu et cb exclu
                        if (malus[l] == minremp && (nbCrdRst = nbCrdRstItr(1, derCrdMin)) > 0)//si il est possible que la ligne saute ?
                        {
                            double probLnSaut = probSaut(calculProbaCrdJr(nbCrdRst, plateau.obtenirNbJ()-1), seuilSaut, puissSaut);
                            int maxLb = 0;
                            int lb = -1;
                            for (int i = 0; i < 4; ++i) if (i != l)
                                {
                                    if (maxLb < derCrd[i] && derCrd[i] < cr)
                                    {
                                        maxLb = derCrd[i];
                                        lb = i;
                                    }
                                }
                            if (lb >= 0)
                            {
                                // nbCrdRstItr(ca, cb) = Nombre de cartes restantes entre ca exclu et cb exclu
                                if (nbCrd[lb] + (Math.Min(1 + (nbCrdRst = nbCrdRstItr(derCrd[lb], cr)), plateau.obtenirNbJ() - 1)) < 6)//On ne rammassera pas !
                                    maluSaut = 0.0f;
                                else
                                {
                                    double probCrdItr = probItrvl(calculProbaCrdJr(nbCrdRst, plateau.obtenirNbJ() - 2), nbCrd[lb]);
                                    maluSaut = malus[l] * probLnSaut * probCrdItr;
                                }
                            }
                            else maluSaut = 0.0f;
                        }
                        else maluSaut = 0.0f;
                        // nbCrdRstItr(ca, cb) = Nombre de cartes restantes entre ca exclu et cb exclu
                        if (nbCrd[l] + (Math.Min(1 + (nbCrdRst = nbCrdRstItr(derCrd[l], cr)), plateau.obtenirNbJ())) < 6)//On ne rammassera pas !
                            malu = 0.0f;
                        else
                        {
                            double probCrdItr = probItrvl(calculProbaCrdJr(nbCrdRst, plateau.obtenirNbJ()-1), nbCrd[l]);
                            malu = malus[l] * probCrdItr;
                        }

                        malu += maluSaut * pondRisq6P;
                    }
                    else malu = minremp * pondPrend;
                    if (malu < min)
                    {
                        res = cr;
                        min = malu;
                    }
                }
                return res;
            }
            return 0;
        }
    }



    class CPioche
    {
        private Random rnd;
        private int nbCartes;
        private byte[] cartes;
        public CPioche()
        {
            rnd = new Random();
            cartes = new byte[104];
            reInitialiser();
        }

        public void reInitialiser()
        {
            nbCartes = 104;
            for (int i = 0; i < cartes.Length; ++i) cartes[i] = 255;
            rnd = new Random(rnd.Next() ^ (new Random().Next()));
        }

        public void remelangerDefausse(IGenJr[] jrs=null)
        {
            for (int i = 0; i < cartes.Length; ++i)
                if(cartes[i]==127)
                {
                    if (jrs!=null) foreach (IGenJr j in jrs) j.piocheMelange((byte)(1+i));
                    cartes[i] = 255;
                    ++nbCartes;
                }
        }

        public int NbCartes() { return nbCartes; }
        public byte piocher()
        {
            if (nbCartes > 0)
            {
                //int num = (int)(Math.Floor(nbCartes * rnd.NextDouble()));
                //if (num >= nbCartes) num = nbCartes - 1;
                int num = rnd.Next(nbCartes);
                for (int i = 0; i < cartes.Length; ++i)
                    if (cartes[i] == 255)
                    {
                        if (num <= 0)
                        {
                            cartes[i] = 0;
                            --nbCartes;
                            return (byte)(1 + i);
                        }
                        else --num;
                    }
            }
            return 0;
        }

        public byte piocher(byte cr)
        {
            if (nbCartes > 0 && 1 <= cr && cr <= 104)
            {
                if(cartes[cr - 1] == 255)
                {
                    cartes[cr - 1] = 0;
                    --nbCartes;
                    return cr;
                }
            }
            return 0;
        }

        public byte piocherAuDessou(byte maxc)
        {
            if (nbCartes > 0)
            {
                maxc = (byte)(maxc - 1);
                int nbCUpper = 0;
                for (int i = 0; i < maxc; ++i)
                    if (cartes[i] > 0) ++nbCUpper;
                if (nbCUpper > 0)
                {
                    int num = rnd.Next(nbCUpper);
                    for (int i = 0; i < maxc; ++i)
                        if (cartes[i] == 255)
                        {
                            if (num <= 0)
                            {
                                cartes[i] = 0;
                                --nbCartes;
                                return (byte)(1 + i);
                            }
                            else --num;
                        }
                }
                return 0;
            }
            return 0;
        }

        static private bool isIn(CInterval[] intrvCs, byte c)
        {
            foreach (CInterval intrvC in intrvCs)
                if (intrvC.a <= c && c <= intrvC.b)
                    return true;
            return false;
        }

        public byte piocher(CInterval[] intrvCs)
        {
            if (nbCartes > 0)
            {
                int nbCIn = 0;
                for (int i = 1; i <= cartes.Length; ++i)
                    if (isIn(intrvCs, (byte)i)) ++nbCIn;

                if (nbCIn > 0)
                {
                    int num = rnd.Next(nbCIn);
                    for (int i = 1; i <= cartes.Length; ++i)
                        if (cartes[i - 1] == 255 && isIn(intrvCs, (byte)i))
                        {
                            if (num <= 0)
                            {
                                cartes[i-1] = 0;
                                --nbCartes;
                                return (byte)(1 + i);
                            }
                            else --num;
                        }
                }
                return 0;
            }
            return 0;
        }

        public void defausser(byte cr)
        {
            if(1<= cr && cr<=104) cartes[cr-1] = 127;
        }
    }
    class CPlateau
    {
        private CPioche pioche;
        private int[] scores;
        private byte[] cardsLine;//4*5

        static public byte nbVache(byte numCrt)
        {
            if (numCrt <= 0 || 104 < numCrt) return 0;

            if (numCrt >= 10)
            {
                if (numCrt == 55) return 7;
                if (numCrt % 10 == 0) return 3;
                if ((numCrt - 5) % 10 == 0) return 2;
                if (numCrt % 10 == numCrt / 10) return 5;
            }
            return 1;
        }

        public CPlateau obtenirFantome()
        {
            CPlateau plat = new CPlateau(scores.Length, null, 66);
            for (int i = 0; i < scores.Length; ++i) plat.scores[i] = scores[i];
            for (int i = 0; i < cardsLine.Length; ++i) plat.cardsLine[i] = cardsLine[i];
            return plat;
        }

        public int obtenirNbJ() { return scores.Length; }

        public CPlateau(int nbJr, CPioche ppioche, int scoreDepart = 66)
        {
            pioche = ppioche;
            scores = new int[nbJr];
            cardsLine = new byte[4 * 5];
            reInitialiser(scoreDepart);
        }

        public void reInitialiser(int scoreDepart = 66)
        {
            for (int i = 0; i < scores.Length; ++i) scores[i] = scoreDepart;
            for (int i = 0; i < cardsLine.Length; ++i)
            {
                if(pioche!=null) pioche.defausser(cardsLine[i]);
                cardsLine[i] = 0;
            }
        }

        public void nouvelleManche(bool melangeEntreManche, IGenJr[] jrs)
        {
            if (pioche != null)
            {
                for (int i = 0; i < cardsLine.Length; ++i)
                {
                    pioche.defausser(cardsLine[i]);
                    cardsLine[i] = 0;
                }
                if (melangeEntreManche) pioche.remelangerDefausse(jrs);
                for (int l = 0; l < 4; ++l)
                {
                    if (pioche.NbCartes() == 0) pioche.remelangerDefausse(jrs);
                    byte cr = pioche.piocher();
                    if (jrs != null) foreach (IGenJr j in jrs) j.jrAction(-1, cr, l);
                    cardsLine[5 * l] = cr;
                }
            }
        }

        public void nouvelleManche(bool melangeEntreManche, byte[] _cardsLine, IGenJr[] jrs)
        {
            if (pioche != null)
            {
                for (int i = 0; i < cardsLine.Length; ++i)
                {
                    pioche.defausser(cardsLine[i]);
                    cardsLine[i] = 0;
                }
                if (melangeEntreManche) pioche.remelangerDefausse(jrs);
                for (int l = 0; l < 4; ++l)
                {
                    if (pioche.NbCartes() == 0) pioche.remelangerDefausse(jrs);
                    byte cr = _cardsLine[5 * l];
                    pioche.piocher(cr);
                    if (jrs != null) foreach (IGenJr j in jrs) j.jrAction(-1, cr, l);
                    cardsLine[5 * l] = cr;
                }
            }
        }

        public int placer(int jr, byte crd, int lig, out int delta)
        {
            int l = ouPlacer(crd);
            if (0 <= l && l < 4)
            {
                if (ajouter(l, crd) == false)
                {
                    delta = malus(l);
                    if (0 <= jr && jr <= scores.Length) scores[jr] -= delta;
                    remplacer(l, crd);
                }
                else delta = 0;
            }
            else if (0 <= lig && lig < 4)
            {
                l = lig;
                delta = malus(l);
                if (0 <= jr && jr <= scores.Length) scores[jr] -= delta;
                remplacer(l, crd);
            }
            else delta = 0;
            return l;
        }

        public int ouPlacer(byte crd)
        {
            int res = -1;
            byte max = 0;

            for(int l=0;l<4;++l)
            {
                int lstCrd = derniereCrd(l);
                if (max < lstCrd && lstCrd < crd)
                {
                    max = (byte)lstCrd;
                    res = l;
                }
            }

            return res;
        }

        public int malus(int lig)
        {
            int res = 0;
            lig *= 5;
            for (int i = 0; i < 5; ++i)res += nbVache(cardsLine[lig + i]);
            return res;
        }

        private void remplacer(int lig, byte crd)
        {
            lig *= 5;
            for (int i = 0; i < 5; ++i)
            {
                if (pioche != null) pioche.defausser(cardsLine[lig + i]);
                cardsLine[lig + i] = 0;
            }
            cardsLine[lig] = crd;
        }

        private bool ajouter(int lig, byte crd)
        {
            lig *= 5;
            int i;
            for (i = 0; i < 5 && (cardsLine[lig + i] > 0); ++i) ;
            if(i<5)
            {
                cardsLine[lig + i] = crd;
                return true;
            }
            else return false;
        }

        /*public byte premiereCrd(int lig)
        {
            return cardsLine[lig*5];
        }*/

        public byte derniereCrd(int lig)
        {
            lig *= 5;
            for(int i = 4; i >= 0; --i)
                if (cardsLine[lig + i]>0)
                    return cardsLine[lig + i];
            return 0;
        }

        public int nbCartes(int lig)
        {
            lig *= 5;
            for (int i = 4; i >= 0; --i)
                if (cardsLine[lig + i] > 0)
                    return (1+i);
            return 0;
        }

        public int Gagnant()
        {
            int res = 0;
            int scoremx = int.MinValue;
            bool finit = false;
            bool matchNull = true;

            for (int i=0;i< scores.Length;++i)
            {
                if(scores[i] < 0) finit = true;
                if (scoremx < scores[i])
                {
                    matchNull = false;
                    scoremx = scores[i];
                    res = i;
                }
                else if (scoremx == scores[i]) matchNull = true;
            }

            return (finit? (matchNull? int.MaxValue:res) : -1);
        }

        //public byte[] GetLines() { return cardsLine; }
        public void CopyLines(byte[] lines) { Array.Copy(cardsLine, lines, cardsLine.Length); }

        public void CalculIntvMin(CInterval[] intervls)
        {
            for (int l = 0; l < 4; ++l)
                intervls[l].a = cardsLine[5 * l] + 1;
            /*{
                if (cardsLine[5 * l + 4] == 0)
                    
                else
                    intervls[l].a = 0;
            }*/
        }

        public void CalculIntvMax(CInterval[] intervls)
        {
            for (int l = 0; l < 4; ++l)
            {
                if (cardsLine[5 * l + 4] == 0)
                {
                    int upp = intervls[l].a;
                    int min = 104 + 1;
                    for (int ol = 0; ol < 4; ++ol)
                        if (ol != l && upp <= intervls[ol].a && intervls[ol].a < min)
                            min = intervls[ol].a;
                    intervls[l].b = min - 1;
                }
                else intervls[l].b = 0;
            }
        }

        public void UpdateIntv(CInterval[] intervls, int lig, byte c)
        {
            c = (byte)(c + 1);
            intervls[lig].a = c;

            int min = 104 + 1;
            if (cardsLine[5 * lig + 4] == 0)
            {
                for (int ol = 0; ol < 4; ++ol)
                    if (ol != lig && c <= intervls[ol].a && intervls[ol].a < min)
                        min = intervls[ol].a;
                intervls[lig].b = min - 1;
            }
            else intervls[lig].b = 0;
        }

        public void RandPlateau(CPioche pioche, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();

            for (int i = 0; i < cardsLine.Length;++i)
                cardsLine[i] = 0;

            CInterval[] intervls = new CInterval[4];
            for (int l = 0; l < 4; ++l)
            {
                byte c = pioche.piocher();
                cardsLine[5 * l] = c;
                intervls[l].a = c;
            }

            int nbC = rnd.Next(0, 16 + 1);

            for (; nbC > 0; --nbC)
            {
                byte c = pioche.piocher();
                if (c > 0)
                {
                    int lig = placer(-1, c, -1, out int delta);
                    UpdateIntv(intervls, lig, c);
                    if (cardsLine[5 * lig + 4] != 0) break;
                }
                else break;
            }

            for (; nbC > 0; --nbC)
            {
                byte c = pioche.piocher(intervls);
                if (c > 0)
                {
                    int lig = placer(-1, c, -1, out int delta);
                    UpdateIntv(intervls, lig, c);
                }
                else break;
            }

            for (; nbC>0; --nbC)
            {
                int nbLigNF = 0;
                byte max = 0;
                for (int l = 0; l < 4; ++l)
                {
                    if (cardsLine[5 * l + 4] == 0)
                    {
                        ++nbLigNF;
                        if (max < cardsLine[5 * l]) max = cardsLine[5 * l];
                    }
                }
                byte c = max > 0 ? pioche.piocherAuDessou(max) : (byte)0;
                if (c != 0 && nbLigNF == 0)
                {
                    int nbl = rnd.Next(nbLigNF);
                    for (int l = 0; l < 4; ++l)
                    {
                        if (cardsLine[5 * l + 4] == 0)
                        {
                            if (nbl > 0) --nbl;
                            else
                            {
                                for (int i = 4; i >= 1; --i)
                                    cardsLine[5 * l + i] = cardsLine[5 * l + i - 1];
                                cardsLine[5 * l] = c;
                            }
                        }
                    }
                }
                else break;
            }
        }
    }
    class CLanceur
    {
        private CPlateau plateau;
        private CPioche pioche;
        private IGenJr[] jrs;
        private byte[] crJr;
        private bool melangeEntreManche;
        private int[] delta;

        public CLanceur(IGenJr[] pjrs, bool pmelangeEntreManche = true)
        {
            jrs = pjrs;
            pioche = new CPioche();
            plateau = new CPlateau(pjrs.Length, pioche);
            delta = new int[pjrs.Length];
            for (int jr = 0; jr < pjrs.Length; ++jr)
            {
                delta[jr] = 0;
                pjrs[jr].setJeu(plateau, pioche, jr);
            }
            crJr = null;
            melangeEntreManche = pmelangeEntreManche;
        }

        public void reInitialiser()
        {
            for (int jr = 0; jr < jrs.Length; ++jr) jrs[jr].reInitialiser();
            plateau.reInitialiser();
            pioche.reInitialiser();
            crJr = null;
        }

        public int realiserPartie()
        {
            int gagnant;
            for (gagnant = -1; gagnant < 0 && realiserManche(); gagnant = plateau.Gagnant()) ;
            return gagnant;
        }

        public bool realiserManche()
        {
            plateau.nouvelleManche(melangeEntreManche, jrs);

            for (int jr = 0; jr < jrs.Length; ++jr)
            {
                for(int i=0;i<10;++i)
                {
                    if (pioche.NbCartes() == 0)pioche.remelangerDefausse(jrs);
                    byte pcrd = pioche.piocher();
                    jrs[jr].recuperer(pcrd);
                }
            }
            bool conti;
            for (conti = true; manche() && conti ; conti = realiserTour()) ;
            return conti;
        }

        public bool realiserTour()
        {
            if (crJr == null || crJr.Length!= jrs.Length)//nouveau tour ???
            {
                crJr = new byte[jrs.Length];
                for (int jr = 0; jr < jrs.Length; ++jr)
                {
                    crJr[jr] = jrs[jr].programmer();
                    if (crJr[jr] == 0)
                    {
                        crJr = null;
                        return false;
                    }
                }
            }
            bool conti = true;
            while(conti)
            {
                int jrAJou = -1;
                byte crdMin = 255;
                for (int jr = 0; jr < jrs.Length; ++jr)
                    if(1 <= crJr[jr] && crJr[jr] <= crdMin)
                    {
                        crdMin = crJr[jr];
                        jrAJou = jr;
                    }
                if (jrAJou < 0)//on a finit le tour
                {
                    crJr = null;
                    conti = false;
                }
                else
                {
                    int lng = plateau.placer(jrAJou, crdMin, -1, out delta[jrAJou]);
                    if (lng < 0)
                    {
                        lng = jrs[jrAJou].choisirLigne(crdMin);
                        if (lng >= 0) lng = plateau.placer(jrAJou, crdMin, lng, out delta[jrAJou]);
                    }

                    if (lng >= 0)
                    {
                        jrs[jrAJou].retirer(crdMin);
                        crJr[jrAJou] = 0;
                        for (int jr = 0; jr < jrs.Length; ++jr)
                            jrs[jr].jrAction(jrAJou, crdMin, lng);
                    }
                    else conti = false;
                }
            }

            for (int jr = 0; jr < jrs.Length; ++jr)
                jrs[jr].scoreDelta(jr, delta);

            return crJr == null;
        }

        private bool manche()
        {
            for (int jr = 0; jr < jrs.Length; ++jr)
                if (jrs[jr].nbCrdMain() <= 0) return false;
            return true;
        }
    }
    class CEtatJeu
    {
        public sbyte[] scores;//suivant nb joueurs
        public byte[] cardsLine;//4*5
        public byte[] cartes;

        public CEtatJeu(int nbJr)
        {
            scores = new sbyte[nbJr]; for(int i=0;i<nbJr;++i) scores[i]=66; cardsLine = new byte[4*5];
            cartes = new byte[104];
            for (int i = 0; i < cartes.Length; ++i) cartes[i] = 255;
        }
        public CEtatJeu(CEtatJeu ej)
        {
            cartes = new byte[104];
            for (int i = 0; i < cartes.Length; ++i) cartes[i] = ej.cartes[i];
            scores = new sbyte[ej.scores.Length];
            for (int i = 0; i < ej.scores.Length; ++i) scores[i] = ej.scores[i];
            cardsLine = new byte[ej.cardsLine.Length];
            for (int i = 0; i < ej.cardsLine.Length; ++i) cardsLine[i] = ej.cardsLine[i];
        }

        static private void carte_mod(byte[] crts, byte[] cj,byte v)
        {
            foreach (byte c in cj)
                if (1 <= c && c <= 104) crts[c - 1] = v;
        }

        public void carte_prog(byte[] cj){carte_mod(cartes, cj,0);}

        private CMMNBCrd ligneMinMaxNbC(int ln, CMMNBCrd mmnbCr = null)
        {
            if (mmnbCr == null) mmnbCr = new CMMNBCrd();
            if (ln < 4)
            {
                ln *= 5;
                mmnbCr.cmin = cardsLine[ln + 0];
                mmnbCr.cmax = 0;
                mmnbCr.nbc = 0;
                mmnbCr.malus = 0;
                for (int cl = 0; cl < 5 && cardsLine[ln + cl]!=0; ++cl)
                {
                    ++mmnbCr.nbc;
                    mmnbCr.cmax = cardsLine[ln + cl];
                    mmnbCr.malus += CTour.nbVache(mmnbCr.cmax);
                }
            }
            else { mmnbCr.malus = mmnbCr.cmin = mmnbCr.cmax = mmnbCr.nbc = 0; mmnbCr.probSaut = 0.0f; }
            return mmnbCr;
        }

        private int nbCrItrvl(int a, int b)
        {
            a--;
            b--;
            int res = 0;

            if (a < 0) a = 0;
            if (b > 103) b = 103;

            for (int i = a; i <= b; ++i)
                if (cartes[i] > 0) ++res;

            return res;
        }

        private int nbVacheItrvl(int a, int b)
        {
            a--;
            b--;
            int res = 0;

            if (a < 0) a = 0;
            if (b > 103) b = 103;

            for (int i = a; i <= b; ++i)
                if (cartes[i] > 0) res += CTour.nbVache((byte)(1+i));

            return res;
        }

        private class CMMNBCrd
        {
            public byte malus;
            public byte cmin;
            public byte cmax;
            public byte nbc;
            public float probSaut;
        }

        private float calcDelta(CMMNBCrd lnMmnbCr, CMMNBCrd[] mmnbCr, byte cr, int d, int nbJ, ref float malPls)
        {
            int lU=-1;
            for (int ln = 0; ln < mmnbCr.Length; ++ln)
                if (mmnbCr[ln] != lnMmnbCr && cr < mmnbCr[ln].cmax && (lU == -1 || mmnbCr[ln].cmax < mmnbCr[lU].cmax))
                    lU = ln;

            int u = nbCrItrvl(lnMmnbCr.cmax + 1, cr - 1);

            int v = 0;
            if (lU>=0) v = nbCrItrvl(cr + 1, mmnbCr[lU].cmax - 1);

            int w = d + u + 1 + v;
            int nbl = 1+(lnMmnbCr.nbc<5 ? (5 - lnMmnbCr.nbc) : 0);

            float probl = (u / ((float)w))*nbJ;

            while(probl>=nbl)
            {
                probl -= nbl;
                nbl = 5;
            }

            malPls = nbVacheItrvl(lnMmnbCr.cmax+1,cr-1)/w;

            return (probl / nbl);
        }

        private int minDistLn(CMMNBCrd[] mmnbCr, int cr, int dMinexclu = 0)
        {
            int minDst = int.MaxValue;
            int sln = -1;
            for(int ln=0;ln< mmnbCr.Length;++ln)
                if(mmnbCr[ln].cmax < cr)
                {
                    int d = (cr - mmnbCr[ln].cmax);
                    if (dMinexclu < d && d < minDst)
                    {
                        minDst = d;
                        sln = ln;
                    }
                }
            return sln;
        }

        private float estimProbMalus(byte cr, CMMNBCrd[] mmnbCr, int d, int nbCrt, int nbJ, float pgs)
        {
            float risq = 0.0f;
            float probPren;
            float malPls=0.0f;

            //float coef = 1.0f - pgs;
            int dMinexclu = 0;
            /*for (int ln = minDistLn(mmnbCr,cr, dMinexclu);ln>=0 && nbJ>0; dMinexclu = (cr-mmnbCr[ln].cmax))
            {
                probPren = calcDelta(mmnbCr[ln], mmnbCr, cr, d, nbJ, ref malPls);
                risq += coef * (probPren * mmnbCr[ln].malus+ malPls);
                --nbJ;
                coef = pgs;
            }*/

            int ln = minDistLn(mmnbCr, cr, dMinexclu);
            if(ln>=0)
            {
                probPren = calcDelta(mmnbCr[ln], mmnbCr, cr, d, nbJ, ref malPls);
                dMinexclu = (cr - mmnbCr[ln].cmax);

                int nln = minDistLn(mmnbCr, cr, dMinexclu);
                if (nln >= 0)
                {
                    risq += (1.0f - pgs) * probPren * (mmnbCr[ln].malus + malPls);
                    byte malus = mmnbCr[ln].malus;
                    byte cmin = mmnbCr[ln].cmin;
                    byte cmax = mmnbCr[ln].cmax;
                    byte nbc = mmnbCr[ln].nbc;
                    float probSaut = mmnbCr[ln].probSaut;

                    for(int c = mmnbCr[ln].cmin-1; c>0 ; --c)
                        if(cartes[c-1]>0)
                        {
                            mmnbCr[ln].cmin = (byte)c;
                            break;
                        }
                    mmnbCr[ln].cmax = mmnbCr[ln].cmin;
                    mmnbCr[ln].nbc = 1;

                    probPren = calcDelta(mmnbCr[nln], mmnbCr, cr, calcNbCrInf(mmnbCr), nbJ-1, ref malPls);
                    mmnbCr[ln].malus = malus;
                    mmnbCr[ln].cmin = cmin;
                    mmnbCr[ln].cmax = cmax;
                    mmnbCr[ln].nbc = nbc;
                    mmnbCr[ln].probSaut = probSaut;
                    risq += pgs * probPren * (mmnbCr[nln].malus + malPls);
                }
                else risq += (probPren * mmnbCr[ln].malus + malPls);
            }
            else
            {
                risq = mmnbCr[0].malus;
                for (ln = 1; ln < mmnbCr.Length; ++ln)
                    if (mmnbCr[ln].malus < risq)
                    {
                        int u = nbCrItrvl(1, cr - 1);
                        //int v = nbCrItrvl(cr + 1, mmnbCr[ln].cmin-1);
                        float p = (u==0?1.0f:(1.0f-(u/((float)d))));
                        risq = mmnbCr[ln].malus * p;
                    }
            }

            return risq;
        }

        private int calcNbCrInf(CMMNBCrd[] mmnbCr)
        {
            byte cmaxmin = byte.MaxValue;
            for (int ln = 0; ln < mmnbCr.Length; ++ln)
                if (mmnbCr[ln].cmax < cmaxmin)
                    cmaxmin = mmnbCr[ln].cmax;

            return nbCrItrvl(1, cmaxmin - 1);
        }

        public string conseuille(byte[] mesCartes)
        {
            string res="";
            float svEsti = float.MaxValue;

            carte_mod(cartes,mesCartes,0);

            CMMNBCrd[] mmnbCr = new CMMNBCrd[5];
            for (int ln = 0; ln < 5; ++ln)
                mmnbCr[ln] = ligneMinMaxNbC(ln);

            int d = calcNbCrInf(mmnbCr);
            int nbCrt = nbCrItrvl(1, 104);
            float pgs = (nbCrt > 0 ? d / ((float)nbCrt) : 0.0f);//prob glob saut

            float cmlprobSaut = 0.0f;
            for (int ln = 0; ln < 5; ++ln)
            {
                float p = (1.1f - mmnbCr[ln].malus / 27.0f);
                mmnbCr[ln].probSaut = pgs * (1.1f - mmnbCr[ln].malus / 27.0f);
                cmlprobSaut += p;
            }
            for (int ln = 0; ln < 5; ++ln)
                mmnbCr[ln].probSaut /= cmlprobSaut;

            for (byte i=0;i<mesCartes.Length && mesCartes[i]>0;++i)
            {
                byte cr = mesCartes[i];
                float estiLoc = estimProbMalus(cr, mmnbCr, d, nbCrt, scores.Length-1, pgs);
                if (estiLoc < svEsti)
                {
                    svEsti = estiLoc;
                    res = cr.ToString();
                }
                else if (estiLoc == svEsti) res += " " + cr.ToString();
            }
            carte_mod(cartes, mesCartes, 255);

            return res;
        }
    }
    class CTour
    {
        public byte[] carteJouees;//suivant nb joueurs
        public CEtatJeu etatFinal;

        static public byte nbVache(byte numCrt)
        {
            if (numCrt <= 0 || 104 < numCrt) return 0;

            if (numCrt >= 10)
            {
                if (numCrt == 55) return 7;
                if (numCrt % 10 == 0) return 3;
                if ((numCrt - 5) % 10 == 0) return 2;
                if (numCrt % 10 == numCrt / 10) return 5;
            }
            return 1;
        }

        private int getIdxLast(int ln)
        {
            if (etatFinal == null || etatFinal.cardsLine == null) return -1;
            int res;
            for (res = 4; (res >= 0 && (etatFinal.cardsLine[ln * 5 + res] == 0)); --res) ;
            return res;
        }

        private int jouerCrtIdx(byte cj)
        {
            if (etatFinal == null || etatFinal.cardsLine == null) return -1;
            int ln, cl, mn;
            ln = cl = mn = -1;
            for (int l = 0; l < 4; ++l)
            {
                int ilst = getIdxLast(l);
                int vl = (ilst < 0 ? 0 : etatFinal.cardsLine[l * 5 + ilst]);
                if (mn < vl && vl < cj)
                {
                    mn = vl;
                    ln = l;
                    cl = ilst;
                }
            }

            return ((mn < 0) ? -1 : (5 * ln + (cl<4 ? (cl+1) : 0)));
        }

        public CTour(byte[] cj, CEtatJeu ej)//On lui passe les cartes qui serront jouées ce tour ci
        {
            carteJouees = new byte[cj.Length];
            for (int i = 0; i < cj.Length; ++i) carteJouees[i] = cj[i];
            etatFinal = new CEtatJeu(ej);
        }

        static public int searchIdxMinCrd(byte[] cj, int imin)
        {
            int idx = -1;
            for (int i = 0; i < cj.Length; ++i)
                if (imin < cj[i] && (idx == -1 || cj[i] < cj[idx]))
                    idx = i;
            return idx;
        }

        public int jouerCartes(byte[] cj, int imin)
        {
            etatFinal.carte_prog(cj);
            int nxt;
            while( (nxt = searchIdxMinCrd(cj, imin))!=-1 )
            {
                byte crt = cj[nxt];
                int idxOu = jouerCrtIdx(crt);
                if (idxOu >= 0)
                {
                    if(etatFinal.cardsLine[idxOu]!=0)
                    {
                        idxOu /= 5;
                        idxOu *= 5;
                        byte malus = 0;
                        for(int i=0;i<5;++i)
                        {
                            malus -= nbVache(etatFinal.cardsLine[idxOu+i]);
                            etatFinal.cardsLine[idxOu + i] = 0;
                        }
                        if (etatFinal.scores != null && nxt < etatFinal.scores.Length)
                            etatFinal.scores[nxt] += (sbyte)malus;
                    }
                    etatFinal.cardsLine[idxOu] = crt;
                    imin = crt;
                }
                else break;
            }
            if (nxt == -1) imin = 127;
            return imin;
        }

        public int jouerCartes(byte[] cj, int imin, int ln)
        {
            int nxt = searchIdxMinCrd(cj, imin);
            if (nxt == -1) return 127;

            byte crt = cj[nxt];
            int idxOu = getIdxLast(ln);
            if (idxOu >= 0)
            {
                if (idxOu >= 4 || crt < etatFinal.cardsLine[5 * ln + idxOu])
                {
                    byte malus = 0;
                    for (int i = 0; i < 5; ++i)
                    {
                        malus += nbVache(etatFinal.cardsLine[5 * ln + i]);
                        etatFinal.cardsLine[5 * ln + i] = 0;
                    }
                    if (etatFinal.scores != null && nxt < etatFinal.scores.Length)
                        etatFinal.scores[nxt] -= (sbyte)malus;
                    idxOu = 0;
                }
                else idxOu += 1;
            }
            else idxOu = 0;
            etatFinal.cardsLine[5*ln + idxOu] = crt;
            imin = crt;

            return imin;
        }
    }

}