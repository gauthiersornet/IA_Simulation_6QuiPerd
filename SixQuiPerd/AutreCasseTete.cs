using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SixQuiPerd
{
    class CSolv
    {
        private byte[] grid;//5l 4c
        private CSolv parent;

        static byte[] cat = { 0, 1, 2, 2, 3, 2, 2, 4, 4, 4, 4 };

        public CSolv() { grid = new byte[4 * 5]; parent = null; }
        public CSolv(CSolv p)
        {
            grid = new byte[4 * 5];
            parent = p;
            for (int i = 0; i < grid.Length; ++i)
                grid[i] = p.grid[i];
        }

        override public bool Equals(Object o)
        {
            if (o is CSolv)
            {
                CSolv csv = (CSolv)o;
                for (int i = 0; i < grid.Length; ++i)
                    if (cat[grid[i]] != cat[csv.grid[i]])//alors test de symétrie
                    {
                        for (int y = 0; y < 5; ++y)
                            for (int x = 0; x < 4; ++x)
                                if (cat[grid[4 * y + (3 - x)]] != cat[csv.grid[4 * y + x]]) return false;
                    }
                return true;
            }
            return false;
        }

        override public int GetHashCode()
        {
            string a = "";
            string b = "";
            for (int y = 0; y < 5; ++y)
                for (int x = 0; x < 4; ++x)
                {
                    a += (char)('A' + cat[grid[4 * y + x]]);
                    b += (char)('A' + cat[grid[4 * y + (3 - x)]]);
                }

            return a.GetHashCode() + b.GetHashCode();
        }

        override public string ToString()
        {
            string r = "";
            for (int y = 0; y < 5; ++y)
            {
                for (int x = 0; x < 4; ++x)
                {
                    byte v = grid[4 * y + x];
                    if (v == 0) r += "_";
                    else r += (char)('A' + v - 1);
                    r += '\t';
                }
                r += "\r\n";
            }
            return r;
        }

        public static string stringAll(CSolv ps)
        {
            int count = 0;
            string r = "";
            while (ps != null)
            {
                ++count;
                r += ps.ToString() + "-\t-\t-\t-distance n°" + count + "\r\n";
                ps = ps.parent;
            }

            return r;
        }

        public static CSolv makeInit()
        {
            CSolv r = new CSolv();
            r.grid[4 * 0 + 0] = 0; r.grid[4 * 0 + 1] = 1; r.grid[4 * 0 + 2] = 1; r.grid[4 * 0 + 3] = 0;
            r.grid[4 * 1 + 0] = 2; r.grid[4 * 1 + 1] = 1; r.grid[4 * 1 + 2] = 1; r.grid[4 * 1 + 3] = 3;
            r.grid[4 * 2 + 0] = 2; r.grid[4 * 2 + 1] = 4; r.grid[4 * 2 + 2] = 4; r.grid[4 * 2 + 3] = 3;
            r.grid[4 * 3 + 0] = 5; r.grid[4 * 3 + 1] = 7; r.grid[4 * 3 + 2] = 8; r.grid[4 * 3 + 3] = 6;
            r.grid[4 * 4 + 0] = 5; r.grid[4 * 4 + 1] = 9; r.grid[4 * 4 + 2] = 10; r.grid[4 * 4 + 3] = 6;
            return r;
        }

        bool isSolved()
        {
            return (grid[4 * 3 + 1] == 1 && grid[4 * 3 + 2] == 1 &&
                    grid[4 * 4 + 1] == 1 && grid[4 * 4 + 2] == 1);

            //return (grid[4 * 0 + 0] == 2);
            //return (grid[4 * 2 + 2] == 8);
        }

        public CSolv move(byte pc, byte dir)
        {
            if (pc == 0) return null;

            CSolv r = new CSolv(this);
            //r.dist++;

            int dx = 0;
            int dy = 0;
            int sy = 0;
            int oy = 1;
            int sx = 0;
            int ox = 1;

            switch (dir)
            {
                case 0:
                    dx = 1;
                    sx = 3;
                    ox = -1;
                    break;
                case 1:
                    dx = -1;
                    break;
                case 2:
                    dy = 1;
                    sy = 4;
                    oy = -1;
                    break;
                case 3:
                    dy = -1;
                    break;
            }

            for (int y = sy; 0 <= y && y < 5; y += oy)
            {
                for (int x = sx; 0 <= x && x < 4; x += ox)
                    if (r.grid[4 * y + x] == pc)
                    {
                        int ny = y + dy;
                        int nx = x + dx;

                        if (0 <= ny && ny < 5 && 0 <= nx && nx < 4 && r.grid[4 * ny + nx] == 0)
                        {
                            r.grid[4 * ny + nx] = pc;
                            r.grid[4 * y + x] = 0;
                        }
                        else return null;
                    }
            }

            return r;
        }

        public static CSolv solve(ref int nbCoupEval)
        {
            LinkedList<CSolv> open = new LinkedList<CSolv>();
            HashSet<CSolv> close = new HashSet<CSolv>();

            {
                CSolv s = makeInit();
                open.AddLast(s);
                close.Add(s);
            }

            CSolv sol = null;
            while (open.First != null && sol == null)
            {
                CSolv s = open.First.Value;
                open.RemoveFirst();


                if (s.isSolved()) sol = s;
                else for (byte pc = 1; pc <= 10; ++pc)
                    {
                        for (byte dir = 0; dir < 4; ++dir)
                        {
                            CSolv dep = s.move(pc, dir);
                            if (dep != null && close.Contains(dep) == false)
                            {
                                open.AddLast(dep);
                                close.Add(dep);
                            }
                        }
                    }
            }

            nbCoupEval = close.Count;
            return sol;
        }
    }


    class CCube3D
    {
        static private bool rsetFrm(byte[] res, byte c, byte ix, byte iy, byte iz, byte pa, byte pb, byte[] frm)
        {
            bool b = true;
            for (int z = -1; z <= 1; ++z)
            {
                int pz = iz + z;
                int fz = 1 + z;
                for (int y = -1; y <= 1; ++y)
                {
                    int py = iy + y;
                    int fy = 1 + y;
                    for (int x = -1; x <= 1; ++x)
                    {
                        int px = ix + x;
                        int fx = 1 + x;
                        if (frm[(fz * 3 + fy) * 3 + fx] == 1)
                        {
                            if ((0 <= pz && pz < c) && (0 <= py && py < c) && (0 <= px && px < c))
                            {
                                int vi = (pz * c + py) * c + px;
                                if (res[vi] == pa) res[vi] = pb;
                                else b = false;
                            }
                            else b = false;
                        }
                    }
                }
            }

            return b;
        }

        static private bool isHole(byte[] res, byte c, byte ix, byte iy, byte iz)
        {
            for (int u = 0; u < 3; ++u)
            {
                for (int xyz = -1; xyz <= 1; xyz += 2)
                {
                    int x = ix + (u == 0 ? xyz : 0);
                    int y = iy + (u == 1 ? xyz : 0);
                    int z = iz + (u == 2 ? xyz : 0);

                    if ((0 <= x && x < c) && (0 <= y && y < c) && (0 <= z && z < c) && res[(z * c + y) * c + x] == 0)
                        return false;
                }
            }
            return true;
        }

        /*static private bool noHoles(byte[] res, byte c, byte ix, byte iy, byte iz)
        {
            for (int z = -2; z <= 0; ++z)if(0 <= (iz+z) && (iz+z) < c)
                for (int y = -2; y <= 2; ++y)if(0 <= (iy+y) && (iy+y) < c)
                    for (int x = -2; x <= 2; ++x)if(0 <= (ix+x) && (ix+x) < c)
                    {
                        if (res[((iz + z) * c + (iy + y)) * c + (ix + x)] == 0 && isHole(res, c, (byte)(ix+x), (byte)(iy+y), (byte)(iz+z)))
                            return false;
                    }
            return true;
        }*/

        static private bool matchForm(byte[] res, byte c, byte ix, byte iy, byte iz, byte[] frm)
        {
            byte b = res[(iz * c + iy) * c + ix];
            if (b == 0) return false;
            for (int z = -1; z <= 1; ++z)
            {
                for (int y = -1; y <= 1; ++y)
                {
                    for (int x = -1; x <= 1; ++x)
                    {
                        if (frm[((1 + z) * 3 + (1 + y)) * 3 + (1 + x)] == 1 &&
                            (((iz + z) < 0 || c <= (iz + z)) || ((iy + y) < 0 || c <= (iy + y)) || ((ix + x) < 0 || c <= (ix + x)) || (res[((iz + z) * c + (iy + y)) * c + (ix + x)] != b))
                            ) return false;
                    }
                }
            }
            return true;
        }

        static private bool noSymetric(byte[] res, byte c, byte ix, byte iy, byte iz, byte[][] formes, int i)
        {
            int dx = (1 + 2 * ix) - c;
            int dy = (1 + 2 * iy) - c;
            int dz = (1 + 2 * iz) - c;

            int uv = i / 2;
            //int u = i / 4;
            //int v = (i / 2)%2;
            int w = (i % 2);
            int frmmir = uv * 2 + (1 - w);
            //int nmFrm;
            //byte[] frm;

            /*byte sc = 3;
            if (uv == 3 || uv == 4) sc = 0;//X
            else if (uv == 0 || uv == 5) sc = 1;//Y
            else if (uv == 1 || uv == 2) sc = 2;//Z*/

            for (int vz = -1; vz <= 1; vz += 2)
                for (int vy = -1; vy <= 1; vy += 2)
                    for (int vx = -1; vx <= 1; vx += 2)
                    {
                        byte px = (byte)((c + dx * vx) / 2);
                        byte py = (byte)((c + dy * vy) / 2);
                        byte pz = (byte)((c + dz * vz) / 2);

                        /*if (sc == 0 && vx < 0) frm = formes[frmmir];
                        else if (sc == 1 && vy < 0) frm = formes[frmmir];
                        else if (sc == 2 && vz < 0) frm = formes[frmmir];
                        else*/ //frm = formes[i];

                        if ((px != ix || py != iy || pz != iz) && (matchForm(res, c, px, py, pz, formes[i]) || matchForm(res, c, px, py, pz, formes[frmmir]))) return false;
                    }

            return true;
        }

        /*static private bool notAHole(byte[] res, byte c, byte ix, byte iy, byte iz, byte[] frm)
        {
        }*/

        /*static private bool noHoles(byte[] res, byte c, byte ix, byte iy, byte iz)
        {

            return true;
        }*/

        /*
         * 0) 13: 64
         * 1) 7 : 96
         * 2) 3 : 48
         * 3) 1 :  8
         * 216 cells
         * 2,1067190168652175306523173339086e+175 combinaisons
        */
        static private bool solv(byte[] res, byte[] sig, byte c, byte ix, byte iy, byte iz, byte p, byte[][] forms, int nbMCb, ArrayList sols, HashSet<CSig> sigs)
        {
            iy += (byte)(ix / c);
            ix = (byte)(ix % c);
            if (iy >= c && sigs.Add(CSig.toSigZYX(sig, c, iz)) == false) return false;
            else sigs = new HashSet<CSig>();
            iz += (byte)(iy / c);
            iy = (byte)(iy % c);

            if (nbMCb == 4 * 0)
            {
                Console.Beep(1000, 100);
                byte[] cpres = new byte[c * c * c];
                res.CopyTo(cpres, 0);
                sols.Add(cpres);
                //CSig.addSet(sigs, sig, c);
                return true;
            }
            else if (iz >= c) return false;
            else if (res[(iz * c + iy) * c + ix] == 0)
            {
                //foreach(byte[] frm in forms)
                for (byte i = 0; i < forms.Length; ++i)
                {
                    /*if (rsetFrm(res, c, ix, iy, iz, 0, p, forms[i]) &&
                        (iz > 0 && (res[(((iz - 1) * c + iy) * c + ix)] == 0 || (iz == (c - 1) && iy > 0 && res[((iz * c + (iy - 1)) * c + ix)] == 0))) == false &&
                        solv(res, c, (byte)(ix + 1), iy, iz, (byte)(p + 1), forms, nbMCb - 4, sols)
                        ) return true;
                    else rsetFrm(res, c, ix, iy, iz, p, 0, forms[i]);*/
                    sig[((iz * c + iy) * c + ix)] = i;
                    if (rsetFrm(res, c, ix, iy, iz, 0, p, forms[i]) &&
                        (iz > 0 && (res[(((iz - 1) * c + iy) * c + ix)] == 0 || (iz == (c - 1) && iy > 0 && res[((iz * c + (iy - 1)) * c + ix)] == 0))) == false)
                        solv(res, sig, c, (byte)(ix + 1), iy, iz, (byte)(p + 1), forms, nbMCb - 4, sols, sigs);
                    rsetFrm(res, c, ix, iy, iz, p, 0, forms[i]);
                    sig[((iz * c + iy) * c + ix)] = (byte)forms.Length;
                }
            }

            return (iz > 0 && (res[(((iz - 1) * c + iy) * c + ix)] == 0 || (iz == (c - 1) && iy > 0 && res[((iz * c + (iy - 1)) * c + ix)] == 0))) == false &&
                    solv(res, sig, c, (byte)(ix + 1), iy, iz, p, forms, nbMCb, sols, sigs);
        }

        static private byte[][] genFroms()
        {
            byte[][] froms = new byte[12][];

            int cf = 1;
            for (int u = 0; u < 3; ++u)
            {
                int cfv = cf;
                for (int v = 0; v <= 1; ++v)
                {
                    cfv = (cfv * 3) % (3 * 3 * 3 - 1);
                    for (int w = -1; w <= 1; w += 2)
                    {
                        byte[] frm = new byte[3 * 3 * 3];
                        for (int xyz = -1; xyz <= 1; ++xyz)
                        {
                            int vi = ((1 * 3 + 1) * 3 + 1) + xyz * cf;
                            frm[vi] = 1;
                        }
                        frm[((1 * 3 + 1) * 3 + 1) + w * cfv] = 1;
                        froms[(u * 2 + v) * 2 + (w < 0 ? 0 : 1)] = frm;
                    }
                }
                cf *= 3;
            }

            return froms;
        }



        static public string cubeToTxt(byte[] cube, byte c)
        {
            string res = "";
            for (int z = 0; z < c; ++z)
            {
                for (int y = 0; y < c; ++y)
                {
                    for (int x = 0; x < c; ++x)
                    {
                        byte cr = cube[(z * c + y) * c + x];
                        if (cr == 0) cr = (byte)'_';
                        else if ((cr - 1) < 26) cr = (byte)('A' + cr - 1);
                        else if ((cr - 1) < 2 * 26) cr = (byte)('a' + cr - 26 - 1);
                        else if ((cr - 1) < (2 * 26 + 10)) cr = (byte)('0' + cr - 2 * (26 - 1));
                        else cr = (byte)'?';
                        res += (char)cr + " ";
                    }
                    res += "\r\n";
                }
                res += "-----------------------------------------------------------------\r\n";
            }
            return res;
        }

        class CSig
        {
            private byte c;
            private byte[] sig;

            public CSig(byte pc)
            {
                c = pc;
                sig = new byte[c * c];
            }

            public override bool Equals(object obj)
            {
                if (obj is CSig)
                {
                    byte[] sigA = ((CSig)obj).sig;
                    if (c > 0 && sigA.Length == c * c && sig.Length == c * c)
                    {
                        byte cm = (byte)(c - 1);

                        bool eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[y * c + x]) eq = false;
                        if (eq) return true;

                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[y * c + (cm - x)]) eq = false;
                        if (eq) return true;

                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[(cm - y) * c + x]) eq = false;
                        if (eq) return true;

                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[(cm - y) * c + (cm - x)]) eq = false;
                        if (eq) return true;



                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[x * c + y]) eq = false;
                        if (eq) return true;

                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[x * c + (cm - y)]) eq = false;
                        if (eq) return true;

                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[(cm - x) * c + y]) eq = false;
                        if (eq) return true;

                        eq = true;
                        for (int y = 0; y < c && eq; ++y)
                            for (int x = 0; x < c && eq; ++x)
                                if (sigA[y * c + x] != sig[(cm - x) * c + (cm - y)]) eq = false;
                        if (eq) return true;
                    }
                }
                return false;
            }

            public override int GetHashCode()
            {
                if (c == 0) return 0;

                string c1, c2, c3, c4, c5, c6, c7, c8;
                c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = "";

                byte cm = (byte)(c - 1);

                for (int y = 0; y < c; ++y)
                    for (int x = 0; x < c; ++x)
                    {
                        c1 += (char)('0' + sig[y * c + x]);
                        c2 += (char)('0' + sig[y * c + (cm - x)]);
                        c3 += (char)('0' + sig[(cm - y) * c + x]);
                        c4 += (char)('0' + sig[(cm - y) * c + (cm - x)]);
                        c5 += (char)('0' + sig[x * c + y]);
                        c6 += (char)('0' + sig[x * c + (cm - y)]);
                        c7 += (char)('0' + sig[(cm - x) * c + y]);
                        c8 += (char)('0' + sig[(cm - x) * c + (cm - y)]);
                    }

                return c1.GetHashCode() + c2.GetHashCode() + c3.GetHashCode() + c4.GetHashCode() + c5.GetHashCode() + c6.GetHashCode() + c7.GetHashCode() + c8.GetHashCode();
            }

            static public CSig toSigZYX(byte[] cub, byte c, int iz)
            {
                if (c == 0 || cub.Length != c * c * c) return null;
                iz *= c * c;
                CSig sig = new CSig(c);
                for (int i = 0; i < c * c; ++i) sig.sig[i] = cub[iz + i];
                return sig;
            }

            static public void addSet(HashSet<CSig> sigs, byte[] sol, byte c)
            {
                if (sol.Length != c * c * c) return;

                CSig sig1, sig2, sig3, sig4, sig5, sig6;
                sig1 = new CSig(c);
                sig2 = new CSig(c);
                sig3 = new CSig(c);
                sig4 = new CSig(c);
                sig5 = new CSig(c);
                sig6 = new CSig(c);
                for (int y = 0; y < c; ++y)
                    for (int x = 0; x < c; ++x)
                    {
                        sig1.sig[y * c + x] = sol[(0 * c + y) * c + x];
                        sig2.sig[y * c + x] = sol[((c - 1) * c + y) * c + x];
                        sig3.sig[y * c + x] = sol[(y * c + 0) * c + x];
                        sig4.sig[y * c + x] = sol[(y * c + (c - 1)) * c + x];
                        sig5.sig[y * c + x] = sol[(y * c + x) * c + 0];
                        sig6.sig[y * c + x] = sol[(y * c + x) * c + (c - 1)];
                    }
                sigs.Add(sig1);
                sigs.Add(sig2);
                sigs.Add(sig3);
                sigs.Add(sig4);
                sigs.Add(sig5);
                sigs.Add(sig6);
            }
        }

        static public ArrayList cubeSolv(byte c = 6)
        {
            ArrayList sols = new ArrayList();
            byte[] res = new byte[c * c * c];
            byte[] sig = new byte[c * c * c];
            byte[][] forms = genFroms();
            for (int i = 0; i < sig.Length; ++i) sig[i] = (byte)forms.Length;
            solv(res, sig, c, 0, 0, 0, 1, forms, c * c * c, sols, new HashSet<CSig>());
            Console.Beep(1000, 1000);
            return sols;
        }
    }
}
