using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SixQuiPerd
{
    public struct SInput
    {
        public byte[] unplayed;//104
        public byte[] hand;//104
        public byte[] lines;//4*5
        static public SInput Zero() { return new SInput { unplayed = new byte[104], hand=new byte[104], lines = new byte[4*5] }; }

        /*public void save(BinaryWriter bw)
        {

        }

        public void load(BinaryReader br)
        {

        }*/
        public void copy(SInput inp)
        {
            Array.Copy(inp.hand, hand, hand.Length);
            Array.Copy(inp.lines, lines, lines.Length);
            Array.Copy(inp.unplayed, unplayed, unplayed.Length);
        }
    }

    public struct SOutput
    {
        public float[] hand;//104
        static public SOutput Zero() { return new SOutput { hand = new float[104] }; }

        /*public void save(BinaryWriter bw)
        {

        }

        public void load(BinaryReader br)
        {

        }*/
    }

    public struct SInpNT
    {
        public float[] unplayed;//104
        public float[] hand;//104
        public float[][] lines;//[104][4*5]

        static public SInpNT Zero()
        {
            float[][] lines = new float[4 * 5][];
            for (int i = 0; i < 4 * 5; ++i)
                lines[i] = new float[104];

            return new SInpNT
            {
                unplayed = new float[104],
                hand = new float[104],
                lines = lines
            };
        }

        static public void RndNT(float[] nt, float min, float max, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();
            for (int i = 0; i < nt.Length; ++i)
                nt[i] = (float)(min + rnd.NextDouble()*(max-min));
        }

        public void RndNT(float min, float max, Random rnd = null)
        {
            if(rnd == null) rnd = new Random();
            RndNT(unplayed, min, max, rnd);
            RndNT(hand, min, max, rnd);
            for (int i = 0; i < 4 * 5; ++i)
                RndNT(lines[i], min, max, rnd);
        }

        public bool Scale(float coeff)
        {
            float sum = unplayed.Sum();
            sum += hand.Sum();
            sum += lines.Sum(l => l.Sum());
            if(-0.0001<sum && sum< 0.0001)
                return false;
            else
            {
                coeff /= sum;
                for (int i = 0; i < unplayed.Length; ++i) unplayed[i] *= coeff;
                for (int i = 0; i < hand.Length; ++i) hand[i] *= coeff;
                for (int i = 0; i < lines.Length; ++i)
                    for (int j = 0; j < lines[i].Length; ++j)
                        lines[i][j] *= coeff;
                return true;
            }
        }

        public void RndNTScale(float min, float max, float coeff, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();
            do {
                RndNT(min, max, rnd);
            } while (Scale(coeff) == false);
        }

        public void save(BinaryWriter bw)
        {
            unplayed.Each(v => bw.Write(v));
            hand.Each(v => bw.Write(v));
            lines.Each(vs => vs.Each(v => bw.Write(v)));
        }

        public SInpNT load(BinaryReader br)
        {
            unplayed = Extentions.Select(104, () => br.ReadSingle()).ToArray();
            hand = Extentions.Select(104, () => br.ReadSingle()).ToArray();
            lines = Extentions.Select(4*5, () => Extentions.Select(104,()=> br.ReadSingle()).ToArray() ).ToArray();
            return this;
        }
    }

    public struct SInpNL
    {
        public SInpNT[] inpNts;
        public SInpNT[] inpNtsN;
        public float[] firstNL;

        static public SInpNL Zero(int cSize)
        {
            SInpNT[] inpNts = new SInpNT[cSize];
            for (int i = 0; i < inpNts.Length; ++i)
                inpNts[i] = SInpNT.Zero();
            SInpNT[] inpNtsN = new SInpNT[cSize];
            for (int i = 0; i < inpNtsN.Length; ++i)
                inpNtsN[i] = SInpNT.Zero();
            float[] firstNL = new float[cSize];
            return new SInpNL { inpNts = inpNts, inpNtsN = inpNtsN, firstNL = firstNL };
        }

        public void RndNT(float min, float max, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();
            for (int i = 0; i < inpNts.Length; ++i)
                inpNts[i].RndNT(min,max,rnd);
            for (int i = 0; i < inpNtsN.Length; ++i)
                inpNtsN[i].RndNT(min, max, rnd);
        }

        public void RndNTScale(float min, float max, float coeff, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();
            for (int i = 0; i < inpNts.Length; ++i)
                inpNts[i].RndNTScale(min, max, coeff, rnd);
            for (int i = 0; i < inpNtsN.Length; ++i)
                inpNtsN[i].RndNTScale(min, max, coeff, rnd);
        }

        /*public void save(BinaryWriter bw)
        {
            unplayed.Each(v => bw.Write(v));
            hand.Each(v => bw.Write(v));
            lines.Each(vs => vs.Each(v => bw.Write(v)));
        }

        public void load(BinaryReader br)
        {
            unplayed = Extentions.Select<float>(104, () => br.ReadSingle()).ToArray();
            hand = Extentions.Select<float>(104, () => br.ReadSingle()).ToArray();
            lines = Extentions.Select<float[]>(4 * 5, () => Extentions.Select<float>(104, () => br.ReadSingle()).ToArray()).ToArray();
        }

        public SInpNT[] inpNts;
        public float[] firstNL;*/

        public void save(BinaryWriter bw)
        {
            bw.Write(inpNts.Length);
            inpNts.Each(v => v.save(bw));
            inpNtsN.Each(v => v.save(bw));
            //firstNL.Each(v => bw.Write(v));
        }

        public SInpNL load(BinaryReader br)
        {
            int ln = br.ReadInt32();
            inpNts = Extentions.Select(ln, () => new SInpNT().load(br)).ToArray();
            inpNtsN = Extentions.Select(ln, () => new SInpNT().load(br)).ToArray();
            firstNL = Extentions.Select(ln, () => 0.0f).ToArray();
            return this;
        }
    }

    public struct SNL
    {
        public float[][] Nts;
        public float[][] NtsN;
        public float[] Ns;

        static public SNL Zero(int cSize, int cSizePrec)
        {
            float[][] Nts = new float[cSize][];
            for (int i = 0; i < Nts.Length; ++i)
                Nts[i] = new float[cSizePrec];
            float[][] NtsN = new float[cSize][];
            for (int i = 0; i < NtsN.Length; ++i)
                NtsN[i] = new float[cSizePrec];
            return new SNL() { Nts = Nts, NtsN = NtsN, Ns = new float[cSize] };
        }

        public void RndNT(float min, float max, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();
            for (int i = 0; i < Nts.Length; ++i)
                for (int j = 0; j < Nts[i].Length; ++j)
                    Nts[i][j] = (float)(min + rnd.NextDouble() * (max - min));
            for (int i = 0; i < NtsN.Length; ++i)
                for (int j = 0; j < NtsN[i].Length; ++j)
                    NtsN[i][j] = (float)(min + rnd.NextDouble() * (max - min));
        }

        public bool Scale(float coeff)
        {
            float sum = Nts.Sum(l => l.Sum());
            if (-0.0001 < sum && sum < 0.0001)
                return false;
            else
            {
                coeff /= sum;
                for (int i = 0; i < Nts.Length; ++i)
                    for (int j = 0; j < Nts[i].Length; ++j)
                        Nts[i][j] *= coeff;
            }

            sum = NtsN.Sum(l => l.Sum());
            if (-0.0001 < sum && sum < 0.0001)
                return false;
            else
            {
                coeff /= sum;
                for (int i = 0; i < NtsN.Length; ++i)
                    for (int j = 0; j < NtsN[i].Length; ++j)
                        NtsN[i][j] *= coeff;
            }

            return true;
        }

        public void RndNTScale(float min, float max, float coeff, Random rnd = null)
        {
            if (rnd == null) rnd = new Random();
            do
            {
                RndNT(min, max, rnd);
            } while (Scale(coeff) == false);
        }

        public void save(BinaryWriter bw)
        {
            bw.Write(Nts.Length);
            bw.Write(Nts.Length>0 ? Nts[0].Length : 0);
            Nts.Each(vs => vs.Each(v => bw.Write(v)));
            NtsN.Each(vs => vs.Each(v => bw.Write(v)));
            //Ns.Each(v => bw.Write(v));
        }

        public SNL load(BinaryReader br)
        {
            int sz = br.ReadInt32();
            int szp = br.ReadInt32();
            Nts = Extentions.Select(sz, () => Extentions.Select(szp, () => br.ReadSingle()).ToArray()).ToArray();
            NtsN = Extentions.Select(sz, () => Extentions.Select(szp, () => br.ReadSingle()).ToArray()).ToArray();
            Ns = Extentions.Select(sz, ()=> 0.0f).ToArray();
            return this;
        }
    }

    class Neural
    {
        protected double age;
        protected double plasticity;
        protected double plastCoef;

        protected SInpNL sInpNL;
        protected SNL[] sNLs;

        public Neural(double _plastCoef, int firstL, int[] hidL, float min, float max, float coeff)
        {
            age = 0.0;
            plasticity = 1.0;
            plastCoef = _plastCoef;

            Random rnd = new Random();
            sInpNL = SInpNL.Zero(firstL);
            sInpNL.RndNTScale(min, max, coeff, rnd);

            int precSize = firstL;
            sNLs = new SNL[hidL.Length+1];
            for (int i = 0; i < hidL.Length; ++i)
            {
                sNLs[i] = SNL.Zero(hidL[i], precSize);
                sNLs[i].RndNTScale(min, max, coeff, rnd);
                precSize = hidL[i];
            }
            sNLs[sNLs.Length - 1] = SNL.Zero(104, precSize);
            sNLs[sNLs.Length - 1].RndNTScale(min, max, coeff, rnd);
        }

        public Neural(string filePath){load(filePath);}

        private void Internal_Inferency(SInput inp)
        {
            for (int i = 0; i < sInpNL.firstNL.Length; ++i)
            {
                SInpNT sInpNT = sInpNL.inpNts[i];
                float sum = 0.0f;
                for (int j = 0; j < inp.unplayed.Length; ++j)
                    if (inp.unplayed[j] > 0) sum += sInpNT.unplayed[j];
                for (int j = 0; j < inp.hand.Length; ++j)
                    if (inp.hand[j] > 0) sum += sInpNT.hand[j];
                for (int j = 0; j < inp.lines.Length; ++j)
                    if (inp.lines[j] > 0) sum += sInpNT.lines[j][inp.lines[j] - 1];
                if (sum < 0.0f) sum = 0.0f;
                else if (sum > 1.0f) sum = 1.0f;
                sInpNL.firstNL[i] = sum;
            }

            float[] precNL = sInpNL.firstNL;
            for (int l = 0; l < sNLs.Length; ++l)
            {
                SNL sNL = sNLs[l];
                for (int j = 0; j < sNL.Ns.Length; ++j)
                {
                    float[] wj = sNL.Nts[j];
                    float sum = 0.0f;
                    for (int k = 0; k < precNL.Length; ++k)
                        sum += wj[k] * precNL[k];
                    if (sum < 0.0f) sum = 0.0f;
                    else if (sum > 1.0f) sum = 1.0f;
                    sNL.Ns[j] = sum;
                }
                precNL = sNL.Ns;
            }
        }

        private void Internal_InferencyRel(SInput inp)
        {
            for (int i = 0; i < sInpNL.firstNL.Length; ++i)
            {
                SInpNT sInpNT = sInpNL.inpNts[i];
                float sum = 0.0f;
                for (int j = 0; j < inp.unplayed.Length; ++j)
                    if (inp.unplayed[j] > 0) sum += sInpNT.unplayed[j];
                    else sum -= sInpNT.unplayed[j];
                for (int j = 0; j < inp.hand.Length; ++j)
                    if (inp.hand[j] > 0) sum += sInpNT.hand[j];
                    else sum -= sInpNT.hand[j];
                for (int j = 0; j < inp.lines.Length; ++j)
                {
                    int kc = inp.lines[j] - 1;
                    for (int k = 0; k < sInpNT.lines[j].Length; ++k)
                        if(k == kc) sum += sInpNT.lines[j][k];
                        else sum -= sInpNT.lines[j][k];
                }
                    /*if (inp.lines[j] > 0) sum += sInpNT.lines[j][inp.lines[j] - 1];
                    else
                    {
                        for(int k=0;k< sInpNT.lines[j].Length;++k)
                            sum -= sInpNT.lines[j][k];
                    }*/
                /*if (sum < 0.0f) sum = 0.0f;
                else if (sum > 1.0f) sum = 1.0f;*/
                if (sum < -1.0f) sum = -1.0f;
                else if (sum > 1.0f) sum = 1.0f;
                sInpNL.firstNL[i] = sum;
            }

            float[] precNL = sInpNL.firstNL;
            for (int l = 0; l < sNLs.Length; ++l)
            {
                SNL sNL = sNLs[l];
                for (int j = 0; j < sNL.Ns.Length; ++j)
                {
                    float[] wj = sNL.Nts[j];
                    float sum = 0.0f;
                    for (int k = 0; k < precNL.Length; ++k)
                        sum += wj[k] * precNL[k];
                    /*if (sum < 0.0f) sum = 0.0f;
                    else if (sum > 1.0f) sum = 1.0f;*/
                    if (sum < -1.0f) sum = -1.0f;
                    else if (sum > 1.0f) sum = 1.0f;
                    sNL.Ns[j] = sum;
                }
                precNL = sNL.Ns;
            }
        }

        private void Internal_InferencySynRel(SInput inp)
        {
            for (int i = 0; i < sInpNL.firstNL.Length; ++i)
            {
                SInpNT sInpNT = sInpNL.inpNts[i];
                SInpNT sInpNTN = sInpNL.inpNtsN[i];
                float sum = 0.0f;
                for (int j = 0; j < inp.unplayed.Length; ++j)
                    if (inp.unplayed[j] > 0) sum += sInpNT.unplayed[j];
                    else sum -= sInpNTN.unplayed[j];
                for (int j = 0; j < inp.hand.Length; ++j)
                    if (inp.hand[j] > 0) sum += sInpNT.hand[j];
                    else sum -= sInpNTN.hand[j];
                for (int j = 0; j < inp.lines.Length; ++j)
                {
                    int kc = inp.lines[j] - 1;
                    for (int k = 0; k < sInpNT.lines[j].Length; ++k)
                        if (k == kc) sum += sInpNT.lines[j][k];
                        else sum -= sInpNTN.lines[j][k];
                }
                /*if (inp.lines[j] > 0) sum += sInpNT.lines[j][inp.lines[j] - 1];
                else
                {
                    for(int k=0;k< sInpNT.lines[j].Length;++k)
                        sum -= sInpNT.lines[j][k];
                }*/
                /*if (sum < 0.0f) sum = 0.0f;
                else if (sum > 1.0f) sum = 1.0f;*/
                if (sum < -1.0f) sum = -1.0f;
                else if (sum > 1.0f) sum = 1.0f;
                sInpNL.firstNL[i] = sum;
            }

            float[] precNL = sInpNL.firstNL;
            for (int l = 0; l < sNLs.Length; ++l)
            {
                SNL sNL = sNLs[l];
                for (int j = 0; j < sNL.Ns.Length; ++j)
                {
                    float[] wj = sNL.Nts[j];
                    float[] wjN = sNL.NtsN[j];
                    float sum = 0.0f;
                    for (int k = 0; k < precNL.Length; ++k)
                        sum += (precNL[k] >=0 ? wj[k] * precNL[k] : wjN[k] * precNL[k]);
                    /*if (sum < 0.0f) sum = 0.0f;
                    else if (sum > 1.0f) sum = 1.0f;*/
                    if (sum < -1.0f) sum = -1.0f;
                    else if (sum > 1.0f) sum = 1.0f;
                    sNL.Ns[j] = sum;
                }
                precNL = sNL.Ns;
            }
        }

        public SOutput Inferency(SInput inp, ref SOutput sout)
        {
            Internal_InferencySynRel(inp);
            Array.Copy(sNLs[sNLs.Length-1].Ns, sout.hand, sout.hand.Length);
            return sout;
        }

        public float BackLearn(SInput inp, SOutput target, float boost = 1.0f)
        {
            float loss = 0.0f;
            boost = Math.Abs(boost);

            float learnCoeff;
            if (plasticity > 0.001)
            {
                learnCoeff = (float)plasticity;
                double newp = plasticity * (1.0 - plastCoef) * boost;
                plasticity -= newp;
            }
            else learnCoeff = 0.001f;
            learnCoeff *= boost;
            age += learnCoeff;

            SNL psNL = sNLs[sNLs.Length - 1];
            for (int i = 0; i < target.hand.Length; ++i)
            {
                psNL.Ns[i] = (target.hand[i] - psNL.Ns[i]) * psNL.Ns[i] * (1.0f - psNL.Ns[i]);
                loss += Math.Abs(psNL.Ns[i]);
            }

            for (int l = sNLs.Length-2; l >= 0; --l)
            {
                SNL sNL = sNLs[l];
                for (int j = 0; j < sNL.Ns.Length; ++j)
                {
                    float delta = 0.0f;
                    for (int k = 0; k < psNL.Ns.Length; ++k)
                    {
                        float w = psNL.Nts[k][j];
                        w += psNL.Ns[k] * sNL.Ns[j] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        psNL.Nts[k][j] = w;
                        delta += psNL.Ns[k] * w;
                    }
                    sNL.Ns[j] = delta * sNL.Ns[j] * (1.0f - sNL.Ns[j]);
                }
                psNL = sNL;
            }

            float[] firstNL = sInpNL.firstNL;
            for (int j = 0; j < firstNL.Length; ++j)
            {
                float delta = 0.0f;
                for (int k = 0; k < psNL.Ns.Length; ++k)
                {
                    float w = psNL.Nts[k][j];
                    w += psNL.Ns[k] * firstNL[j] * learnCoeff;
                    if (w < -1.0f) w = -1.0f;
                    else if (w > 1.0f) w = 1.0f;
                    psNL.Nts[k][j] = w;
                    delta += psNL.Ns[k] * w;
                }
                firstNL[j] = delta * firstNL[j] * (1.0f - firstNL[j]);
            }

            for (int i = 0; i < sInpNL.firstNL.Length; ++i)
            {
                SInpNT sInpNT = sInpNL.inpNts[i];
                for (int j = 0; j < inp.unplayed.Length; ++j)
                    if (inp.unplayed[j] > 0)
                    {
                        //if(inp.unplayed[j]>0)sum += sInpNT.unplayed[j];
                        float w = sInpNT.unplayed[j];
                        w += firstNL[i] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNT.unplayed[j] = w;
                    }
                for (int j = 0; j < inp.hand.Length; ++j)
                    if (inp.hand[j] > 0)
                    {
                        //if (inp.hand[j] > 0) sum += sInpNT.hand[j];
                        float w = sInpNT.hand[j];
                        w += sInpNT.unplayed[j] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNT.hand[j] = w;
                    }
                for (int j = 0; j < inp.lines.Length; ++j)
                    if (inp.lines[j] > 0)
                    {
                        //if (inp.lines[j] > 0) sum += sInpNT.lines[j][inp.lines[j] - 1];
                        float w = sInpNT.lines[j][inp.lines[j] - 1];
                        w += firstNL[i] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNT.lines[j][inp.lines[j] - 1] = w;
                    }
            }

            return loss;
        }

        public float BackLearnSynRel(SInput inp, SOutput target, float boost = 1.0f)
        {
            float loss = 0.0f;
            boost = Math.Abs(boost);

            float learnCoeff;
            if (plasticity > 0.001)
            {
                learnCoeff = (float)plasticity;
                double newp = plasticity * (1.0 - plastCoef) * boost;
                plasticity -= newp;
            }
            else learnCoeff = 0.001f;
            learnCoeff *= boost;
            age += learnCoeff;

            SNL psNL = sNLs[sNLs.Length - 1];
            for (int i = 0; i < target.hand.Length; ++i)
            {
                psNL.Ns[i] = (target.hand[i] - psNL.Ns[i]) * psNL.Ns[i] * (1.0f - psNL.Ns[i]);
                loss += Math.Abs(psNL.Ns[i]);
            }

            for (int l = sNLs.Length - 2; l >= 0; --l)
            {
                SNL sNL = sNLs[l];
                for (int j = 0; j < sNL.Ns.Length; ++j)
                {
                    float delta = 0.0f;
                    for (int k = 0; k < psNL.Ns.Length; ++k)
                        if(psNL.Ns[k] >= 0.0f)
                        {
                            float w = psNL.Nts[k][j];
                            w += psNL.Ns[k] * sNL.Ns[j] * learnCoeff;
                            if (w < -1.0f) w = -1.0f;
                            else if (w > 1.0f) w = 1.0f;
                            psNL.Nts[k][j] = w;
                            delta += psNL.Ns[k] * w;
                        }
                        else
                        {
                            float w = psNL.NtsN[k][j];
                            w += psNL.Ns[k] * sNL.Ns[j] * learnCoeff;
                            if (w < -1.0f) w = -1.0f;
                            else if (w > 1.0f) w = 1.0f;
                            psNL.NtsN[k][j] = w;
                            delta += psNL.Ns[k] * w;
                        }
                    sNL.Ns[j] = delta * sNL.Ns[j] * (1.0f - sNL.Ns[j]);
                }
                psNL = sNL;
            }

            float[] firstNL = sInpNL.firstNL;
            for (int j = 0; j < firstNL.Length; ++j)
            {
                float delta = 0.0f;
                for (int k = 0; k < psNL.Ns.Length; ++k)
                    if(psNL.Ns[k] >= 0.0f)
                    {
                        float w = psNL.Nts[k][j];
                        w += psNL.Ns[k] * firstNL[j] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        psNL.Nts[k][j] = w;
                        delta += psNL.Ns[k] * w;
                    }
                    else
                    {
                        float w = psNL.NtsN[k][j];
                        w += psNL.Ns[k] * firstNL[j] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        psNL.NtsN[k][j] = w;
                        delta += psNL.Ns[k] * w;
                    }
                firstNL[j] = delta * firstNL[j] * (1.0f - firstNL[j]);
            }

            for (int i = 0; i < sInpNL.firstNL.Length; ++i)
            {
                SInpNT sInpNT = sInpNL.inpNts[i];
                SInpNT sInpNTN = sInpNL.inpNtsN[i];
                for (int j = 0; j < inp.unplayed.Length; ++j)
                    if (inp.unplayed[j] > 0)
                    {
                        //if(inp.unplayed[j]>0)sum += sInpNT.unplayed[j];
                        float w = sInpNT.unplayed[j];
                        w += firstNL[i] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNT.unplayed[j] = w;
                    }
                    else
                    {
                        //if(inp.unplayed[j]>0)sum += sInpNT.unplayed[j];
                        float w = sInpNTN.unplayed[j];
                        w -= firstNL[i] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNTN.unplayed[j] = w;
                    }

                for (int j = 0; j < inp.hand.Length; ++j)
                    if (inp.hand[j] > 0)
                    {
                        //if (inp.hand[j] > 0) sum += sInpNT.hand[j];
                        float w = sInpNT.hand[j];
                        w += sInpNT.unplayed[j] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNT.hand[j] = w;
                    }
                    else
                    {
                        //if (inp.hand[j] > 0) sum += sInpNT.hand[j];
                        float w = sInpNTN.hand[j];
                        w -= sInpNT.unplayed[j] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNTN.hand[j] = w;
                    }

                for (int j = 0; j < inp.lines.Length; ++j)
                    if (inp.lines[j] > 0)
                    {
                        //if (inp.lines[j] > 0) sum += sInpNT.lines[j][inp.lines[j] - 1];
                        float w = sInpNT.lines[j][inp.lines[j] - 1];
                        w += firstNL[i] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNT.lines[j][inp.lines[j] - 1] = w;
                    }
                    else
                    {
                        //if (inp.lines[j] > 0) sum += sInpNT.lines[j][inp.lines[j] - 1];
                        float w = sInpNTN.lines[j][inp.lines[j] - 1];
                        w -= firstNL[i] * learnCoeff;
                        if (w < -1.0f) w = -1.0f;
                        else if (w > 1.0f) w = 1.0f;
                        sInpNTN.lines[j][inp.lines[j] - 1] = w;
                    }
            }

            return loss;
        }

        public float Trainning(IGenJr tutor, int nbJr, int nbTurn)
        {
            float loss = 0.0f;
            SInput inp = SInput.Zero();
            SOutput sOut = SOutput.Zero();

            CPioche pioche = new CPioche();
            CPlateau plateau = new CPlateau(nbJr, pioche);

            Random rnd = new Random();
            tutor.setJeu(plateau, pioche, nbJr);

            for (; nbTurn>0;--nbTurn)
            {
                tutor.reInitialiser();
                rnd = new Random(rnd.Next() ^ (new Random()).Next());

                pioche.reInitialiser();
                for (int i = 0; i < inp.hand.Length; ++i)
                {
                    inp.hand[i] = 0;
                    sOut.hand[i] = -1.0f;
                }

                plateau.RandPlateau(pioche, rnd);
                for (int i = 0; i < 10; ++i)
                {
                    byte c = pioche.piocher();
                    inp.hand[c - 1] = 255;
                    tutor.recuperer(c);
                }
                plateau.CopyLines(inp.lines);

                byte chx = tutor.programmer();
                sOut.hand[chx] = 1.0f;

                Internal_Inferency(inp);
                loss = BackLearnSynRel(inp, sOut, 1.0f);
            }

            return loss;
        }

        public void save(string filePath)
        {
            using (BinaryWriter bw = new BinaryWriter(File.OpenWrite(filePath)))
            {
                bw.Write(age);
                bw.Write(plasticity);
                bw.Write(plastCoef);

                sInpNL.save(bw);

                bw.Write(sNLs.Length);
                sNLs.Each(snl => snl.save(bw));
                bw.Close();
            }
        }

        public void load(string filePath)
        {
            using (BinaryReader br = new BinaryReader(File.OpenRead(filePath)))
            {
                age = br.ReadDouble();
                plasticity = br.ReadDouble();
                plastCoef = br.ReadDouble();

                sInpNL.load(br);

                int len = br.ReadInt32();
                sNLs = Extentions.Select(len, () => new SNL().load(br)).ToArray();
            }
        }
    }

    class BasicDeepLIA : CProgCptJRunSup, IGenJr
    {
        protected Neural _nn;

        int nbHC;
        SInput inp;
        SOutput outp;
        SOutput target;
        byte chxCrd;

        float seuilDeci = 0.05f;

        public long Errors { get; private set; }

        public BasicDeepLIA(double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
            : base(ppondPrend, ppondRisq6P)
        {
            _nn = new Neural(0.999, 3 * 104, new int[] { 2 * 104, (3 * 104) / 2 }, -0.5f, 0.5f, 0.7f);
            plateau = null;
            jrNbr = 0;
            rnd = new Random();
            inp = SInput.Zero();
            nbHC = 0;
            for (int i = 0; i < 104; ++i) inp.hand[i] = 0;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;
            outp = SOutput.Zero();
            target = SOutput.Zero();
            Errors = 0;
        }

        public BasicDeepLIA(string filePath, double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
            : base(ppondPrend, ppondRisq6P)
        {
            _nn = new Neural(filePath);
            plateau = null;
            jrNbr = 0;
            rnd = new Random();
            inp = SInput.Zero();
            nbHC = 0;
            for (int i = 0; i < 104; ++i) inp.hand[i] = 0;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;
            outp = SOutput.Zero();
            target = SOutput.Zero();
            Errors = 0;
        }

        public void SaveNN(string filePath)
        {
            _nn.save(filePath);
        }

        override public void jrAction(int jr, byte cr, int ln)
        {
            if (1 <= cr && cr <= 104) inp.unplayed[cr - 1] = 0;
            base.jrAction(jr, cr, ln);
        }

        override public int nbCrdMain() { return nbHC; }

        override public void piocheMelange()
        {
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;
            base.piocheMelange();
        }

        override public void piocheMelange(byte cr)
        {
            if (1 <= cr && cr <= 104) inp.unplayed[cr - 1] = 255;
            base.piocheMelange(cr);
        }

        override public byte programmer()
        {
            if (nbHC > 0)
            {
                plateau.CopyLines(inp.lines);
                SOutput sout = _nn.Inferency(inp, ref outp);

                float max = 0.0f;
                for (int i = 0; i < inp.hand.Length; ++i)
                    if (inp.hand[i] > 0)
                    {
                        if (sout.hand[i] > max)
                        {
                            max = sout.hand[i];
                            chxCrd = (byte)i;
                        }
                    }

                if (max <= seuilDeci)
                {
                    /*//Random rnd = new Random(DateTime.Now.Millisecond);
                    int num = (int)(Math.Floor(nbHC * rnd.NextDouble()));
                    if (num >= nbHC) num = nbHC - 1;
                    for (int i = 0; i < inp.hand.Length; ++i)
                        if (inp.hand[i] > 0)
                        {
                            chxCrd = (byte)i;
                            if (num == 0)break;
                            else --num;
                        }*/
                    chxCrd = 0;
                    return base.programmer();
                }

                return (byte)(chxCrd + 1);
            }
            return 0;
        }

        override public void recuperer(byte c)
        {
            if (1 <= c && c <= 104 && inp.hand[c - 1] == 0)
            {
                ++nbHC;
                inp.unplayed[c - 1] = 0;
                inp.hand[c - 1] = 255;
            }
            base.recuperer(c);
        }

        override public void reInitialiser()
        {
            inp = SInput.Zero();
            for (int i = 0; i < 104; ++i) inp.hand[i] = 0;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;

            nbHC = 0;
            base.reInitialiser();
        }

        override public void retirer(byte c)
        {
            if (nbHC > 0 && 1 <= c && c <= 104 && inp.hand[c - 1] > 0)
            {
                --nbHC;
                inp.hand[c - 1] = 0;
            }
            base.retirer(c);
        }

        private readonly float MoyMalus = 8.6538461538461538461538461538461f;
        private readonly float LostBoost = -1.5f;
        private readonly float WinBoost = 2.0f;
        private readonly float Boost = 0.5f;

        override public void scoreDelta(int nmJr, int[] delta)
        {
            if (chxCrd > 0)
            {
                int ownDelta = delta[nmJr];
                float boost;

                if (ownDelta > 0)
                {
                    double avg = delta.Average();
                    if (ownDelta > avg)
                    {
                        if (ownDelta > MoyMalus) boost = LostBoost;
                        else
                        {
                            boost = (LostBoost * ownDelta) / MoyMalus;
                            if (boost > -1.0f) boost = -1.0f;
                        }
                    }
                    else boost = (float)((LostBoost * ownDelta) / avg);
                }
                else
                {
                    int max = delta.Max();
                    if (max > 0)
                    {
                        if (max > MoyMalus) boost = WinBoost;
                        else boost = (WinBoost * max) / MoyMalus;
                    }
                    else boost = Boost;
                }

                int idxChxCr = chxCrd - 1;
                if (boost > 0.0f)
                {
                    Errors -= 2;
                    for (int i = 0; i < target.hand.Length; ++i)
                        target.hand[i] = (i == idxChxCr ? 1.0f : -1.0f);
                }
                else if (boost < 0.0f)
                {
                    Errors += 3;
                    for (int i = 0; i < target.hand.Length; ++i)
                        target.hand[i] = (i == idxChxCr ? -1.0f : 1.0f);
                }

                _nn.BackLearnSynRel(inp, target, boost);
            }
            else Errors += 1;
        }
    }

    class QDeepLIA : CProgCptJRunSup, IGenJr
    {
        protected Neural _nn;

        int nbHC;
        SInput inpSv;
        SInput inp;
        bool changed;
        SOutput outp;
        SOutput target;
        byte chxCrd;

        float seuilDeci = 0.05f;

        float alpha = 0.01f;
        float rabais = 0.95f;

        public long Errors { get; private set; }

        public QDeepLIA(double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
            : base(ppondPrend, ppondRisq6P)
        {
            _nn = new Neural(0.999, 3 * 104, new int[] { 2 * 104, (3 * 104) / 2 }, -0.5f, 0.5f, 0.7f);
            plateau = null;
            jrNbr = 0;
            rnd = new Random();
            inpSv = SInput.Zero();
            inp = SInput.Zero();
            nbHC = 0;
            for (int i = 0; i < 104; ++i) inp.hand[i] = 0;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;
            outp = SOutput.Zero();
            target = SOutput.Zero();
            changed = true;
            Errors = 0;
        }

        public QDeepLIA(string filePath, double ppondPrend = 2.5f, double ppondRisq6P = 0.15f)
            : base(ppondPrend, ppondRisq6P)
        {
            _nn = new Neural(filePath);
            plateau = null;
            jrNbr = 0;
            rnd = new Random();
            inp = SInput.Zero();
            nbHC = 0;
            for (int i = 0; i < 104; ++i) inp.hand[i] = 0;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;
            outp = SOutput.Zero();
            target = SOutput.Zero();
            changed = true;
            Errors = 0;
        }

        public void SaveNN(string filePath)
        {
            _nn.save(filePath);
        }

        override public void jrAction(int jr, byte cr, int ln)
        {
            if (1 <= cr && cr <= 104)
            {
                inp.unplayed[cr - 1] = 0;
                changed = true;
            }
            base.jrAction(jr, cr, ln);
        }

        override public int nbCrdMain() { return nbHC; }

        override public void piocheMelange()
        {
            changed = true;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;
            base.piocheMelange();
        }

        override public void piocheMelange(byte cr)
        {
            if (1 <= cr && cr <= 104)
            {
                changed = true;
                inp.unplayed[cr - 1] = 255;
            }
            base.piocheMelange(cr);
        }

        override public byte programmer()
        {
            if (nbHC > 0)
            {
                plateau.CopyLines(inp.lines);
                inpSv.copy(inp);

                if (changed)
                {
                    _nn.Inferency(inp, ref outp);
                    changed = false;
                }

                float max = 0.0f;
                for (int i = 0; i < inp.hand.Length; ++i)
                    if (inp.hand[i] > 0)
                    {
                        if (outp.hand[i] > max)
                        {
                            max = outp.hand[i];
                            chxCrd = (byte)i;
                        }
                    }

                if (max <= seuilDeci)
                {
                    /*//Random rnd = new Random(DateTime.Now.Millisecond);
                    int num = (int)(Math.Floor(nbHC * rnd.NextDouble()));
                    if (num >= nbHC) num = nbHC - 1;
                    for (int i = 0; i < inp.hand.Length; ++i)
                        if (inp.hand[i] > 0)
                        {
                            chxCrd = (byte)i;
                            if (num == 0)break;
                            else --num;
                        }*/
                    chxCrd = 0;
                    return base.programmer();
                }

                return (byte)(chxCrd + 1);
            }
            return 0;
        }

        override public void recuperer(byte c)
        {
            if (1 <= c && c <= 104 && inp.hand[c - 1] == 0)
            {
                ++nbHC;
                inp.unplayed[c - 1] = 0;
                inp.hand[c - 1] = 255;
            }
            base.recuperer(c);
        }

        override public void reInitialiser()
        {
            inp = SInput.Zero();
            for (int i = 0; i < 104; ++i) inp.hand[i] = 0;
            for (int i = 0; i < 104; ++i) inp.unplayed[i] = 255;

            nbHC = 0;
            base.reInitialiser();
        }

        override public void retirer(byte c)
        {
            if (nbHC > 0 && 1 <= c && c <= 104 && inp.hand[c - 1] > 0)
            {
                --nbHC;
                inp.hand[c - 1] = 0;
            }
            base.retirer(c);
        }

        private readonly float MoyMalus = 8.6538461538461538461538461538461f;
        private readonly float LostBoost = -1.0f;
        private readonly float WinBoost = 1.0f;
        private readonly float Boost = 0.5f;

        override public void scoreDelta(int nmJr, int[] delta)
        {
            if (chxCrd > 0)
            {
                int ownDelta = delta[nmJr];
                /*float boost;

                if (ownDelta > 0)
                {
                    double avg = delta.Average();
                    if (ownDelta > avg)
                    {
                        if (ownDelta > MoyMalus) boost = LostBoost;
                        else
                        {
                            boost = (LostBoost * ownDelta) / MoyMalus;
                            if (boost > -1.0f) boost = -1.0f;
                        }
                    }
                    else boost = (float)((LostBoost * ownDelta) / avg);
                }
                else
                {
                    int max = delta.Max();
                    if (max > 0)
                    {
                        if (max > MoyMalus) boost = WinBoost;
                        else boost = (WinBoost * max) / MoyMalus;
                    }
                    else boost = Boost;
                }*/

                /*int idxChxCr = chxCrd - 1;
                if (boost > 0.0f)
                {
                    Errors -= 2;
                    for (int i = 0; i < target.hand.Length; ++i)
                        target.hand[i] = (i == idxChxCr ? 1.0f : -1.0f);
                }
                else if (boost < 0.0f)
                {
                    Errors += 3;
                    for (int i = 0; i < target.hand.Length; ++i)
                        target.hand[i] = (i == idxChxCr ? -1.0f : 1.0f);
                }*/

                float recomp;
                if (ownDelta > 0)
                {
                    double avg = delta.Average();
                    if (ownDelta > avg)
                    {
                        if (ownDelta > MoyMalus) recomp = LostBoost;
                        else
                        {
                            recomp = (LostBoost * ownDelta) / MoyMalus;
                            if (recomp > -1.0f) recomp = -1.0f;
                        }
                    }
                    else recomp = (float)((LostBoost * ownDelta) / avg);
                }
                else
                {
                    int max = delta.Max();
                    if (max > 0)
                    {
                        if (max > MoyMalus) recomp = WinBoost;
                        else recomp = (WinBoost * max) / MoyMalus;
                    }
                    else recomp = Boost;
                }

                //Ajouter à la récompense le fait d'avoir gagné ou pas...
                //...

                plateau.CopyLines(inp.lines);
                _nn.Inferency(inp, ref outp);
                changed = false;

                float fmax = 0.0f;
                for (int i = 0; i < outp.hand.Length; ++i)
                    if (inp.hand[i] > 0.0f && outp.hand[i] > fmax)
                        fmax = outp.hand[i];

                _nn.Inferency(inpSv, ref target);

                if(0<=chxCrd && chxCrd<104)
                {
                    target.hand[chxCrd] += alpha * (recomp + rabais * fmax - target.hand[chxCrd]);
                    float mx = target.hand.Max();
                    float mn = target.hand.Min();
                    if (mn != mx)
                    {
                        float dltMnMx = mx - mn;
                        for (int i = 0; i < target.hand.Length; ++i)
                            target.hand[i] = mn + ((target.hand[i] - mn) / dltMnMx);
                    }
                    _nn.BackLearnSynRel(inp, target);
                }
            }
            else Errors += 1;
        }
    }
}
