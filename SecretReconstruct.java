import java.nio.file.*;
import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.regex.*;

public class SecretReconstruct {

    static class Share {
        int x;
        int base;
        String valueStr;
        BigInteger y;
    }
    static class InputData {
        int n, k;
        List<Share> shares = new ArrayList<>();
    }

    static InputData parseJsonTight(String json) {
        InputData data = new InputData();

        Pattern pNK = Pattern.compile("\"keys\"\\s*:\\s*\\{[^}]*\"n\"\\s*:\\s*(\\d+)\\s*,\\s*\"k\"\\s*:\\s*(\\d+)[^}]*\\}");
        Matcher mNK = pNK.matcher(json);
        if (!mNK.find()) throw new RuntimeException("Could not parse keys.n and keys.k");
        data.n = Integer.parseInt(mNK.group(1));
        data.k = Integer.parseInt(mNK.group(2));

        Pattern pShare = Pattern.compile("\"(\\d+)\"\\s*:\\s*\\{\\s*\"base\"\\s*:\\s*\"(\\d+)\"\\s*,\\s*\"value\"\\s*:\\s*\"([0-9A-Za-z]+)\"\\s*\\}");
        Matcher ms = pShare.matcher(json);
        while (ms.find()) {
            Share s = new Share();
            s.x = Integer.parseInt(ms.group(1));
            s.base = Integer.parseInt(ms.group(2));
            s.valueStr = ms.group(3);
            data.shares.add(s);
        }
        if (data.shares.size() != data.n) {
            throw new RuntimeException("n does not match number of parsed shares (" + data.n + " vs " + data.shares.size() + ")");
        }
        return data;
    }

    static BigInteger parseInBase(String digits, int base) {
        BigInteger B = BigInteger.valueOf(base);
        BigInteger acc = BigInteger.ZERO;
        for (int i = 0; i < digits.length(); i++) {
            char c = digits.charAt(i);
            int v;
            if (c >= '0' && c <= '9') v = c - '0';
            else if (c >= 'a' && c <= 'z') v = 10 + (c - 'a');
            else if (c >= 'A' && c <= 'Z') v = 10 + (c - 'A');
            else throw new IllegalArgumentException("Invalid digit: " + c);
            if (v >= base) throw new IllegalArgumentException("Digit out of range for base " + base + ": " + c);
            acc = acc.multiply(B).add(BigInteger.valueOf(v));
        }
        return acc;
    }

    static BigInteger mod(BigInteger a, BigInteger p) {
        a = a.mod(p);
        return a.signum() < 0 ? a.add(p) : a;
    }
    static BigInteger modAdd(BigInteger a, BigInteger b, BigInteger p) { return mod(a.add(b), p); }
    static BigInteger modMul(BigInteger a, BigInteger b, BigInteger p) { return mod(a.multiply(b), p); }
    static BigInteger modInv(BigInteger a, BigInteger p) { return a.modInverse(p); }
    static BigInteger choosePrimeAbove(BigInteger bound) { return bound.nextProbablePrime(); }

    // f(0) via Lagrange-at-zero with the k points in xs, ys
    static BigInteger secretAtZeroFromSubset(int[] xs, BigInteger[] ys, BigInteger p) {
        int k = xs.length;
        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < k; i++) {
            BigInteger numer = BigInteger.ONE;
            BigInteger denom = BigInteger.ONE;
            BigInteger xi = BigInteger.valueOf(xs[i]);
            for (int j = 0; j < k; j++) {
                if (i == j) continue;
                BigInteger xj = BigInteger.valueOf(xs[j]);
                numer = modMul(numer, mod(xj.negate(), p), p);     
                denom = modMul(denom, mod(xi.subtract(xj), p), p);  
            }
            BigInteger li0 = modMul(numer, modInv(denom, p), p);
            sum = modAdd(sum, modMul(ys[i], li0, p), p);
        }
        return sum;
    }

   
    static BigInteger evalAtXFromSubset(int[] xs, BigInteger[] ys, int X, BigInteger p) {
        int k = xs.length;
        BigInteger x = BigInteger.valueOf(X);
        BigInteger sum = BigInteger.ZERO;
        for (int i = 0; i < k; i++) {
            BigInteger numer = BigInteger.ONE;
            BigInteger denom = BigInteger.ONE;
            BigInteger xi = BigInteger.valueOf(xs[i]);
            for (int j = 0; j < k; j++) {
                if (i == j) continue;
                BigInteger xj = BigInteger.valueOf(xs[j]);
                numer = modMul(numer, mod(x.subtract(xj), p), p);   
                denom = modMul(denom, mod(xi.subtract(xj), p), p);   
            }
            BigInteger Li = modMul(numer, modInv(denom, p), p);
            sum = modAdd(sum, modMul(ys[i], Li, p), p);
        }
        return sum;
    }

   
    static void combinations(int n, int k, IntConsumer action) {
        int[] idx = new int[k];
        for (int i = 0; i < k; i++) idx[i] = i;
        while (true) {
            action.accept(Arrays.copyOf(idx, k));
            int i = k - 1;
            while (i >= 0 && idx[i] == i + n - k) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
        }
    }
    interface IntConsumer { void accept(int[] comb); }

   
    static class Best {
        int agree = -1;
        BigInteger secret = null;
        int[] comb = null;
        final List<int[]> ties = new ArrayList<>();
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            System.err.println("Usage: java SecretReconstruct <input.json>");
            System.exit(1);
        }
        String json = Files.readString(Paths.get(args[0]));

        InputData data = parseJsonTight(json);
        data.shares.sort(Comparator.comparingInt(s -> s.x));

        BigInteger maxVal = BigInteger.ZERO;
        for (Share s : data.shares) {
            s.y = parseInBase(s.valueStr, s.base);
            if (s.y.compareTo(maxVal) > 0) maxVal = s.y;
            if (BigInteger.valueOf(s.x).compareTo(maxVal) > 0) maxVal = BigInteger.valueOf(s.x);
        }

        BigInteger p = choosePrimeAbove(maxVal.add(BigInteger.valueOf(5)));

        int n = data.n, k = data.k;
        List<Integer> xsAll = new ArrayList<>();
        List<BigInteger> ysAll = new ArrayList<>();
        for (Share s : data.shares) {
            xsAll.add(s.x);
            ysAll.add(mod(s.y, p));
        }

        final Best best = new Best();

        combinations(n, k, comb -> {
            int[] xs = new int[k];
            BigInteger[] ys = new BigInteger[k];
            for (int i = 0; i < k; i++) {
                xs[i] = xsAll.get(comb[i]);
                ys[i] = ysAll.get(comb[i]);
            }
            BigInteger s0 = secretAtZeroFromSubset(xs, ys, p);

            int agree = 0;
            for (int t = 0; t < n; t++) {
                int xT = xsAll.get(t);
                BigInteger yT = ysAll.get(t);
                BigInteger yHat = evalAtXFromSubset(xs, ys, xT, p);
                if (yHat.equals(yT)) agree++;
            }

            if (agree > best.agree) {
                best.agree = agree;
                best.secret = s0;
                best.comb = Arrays.copyOf(comb, comb.length);
                best.ties.clear();
                best.ties.add(best.comb);
            } else if (agree == best.agree) {
                best.ties.add(Arrays.copyOf(comb, comb.length));
            }
        });

        int[] xsBest = new int[k];
        BigInteger[] ysBest = new BigInteger[k];
        for (int i = 0; i < k; i++) {
            xsBest[i] = xsAll.get(best.comb[i]);
            ysBest[i] = ysAll.get(best.comb[i]);
        }

        List<String> mismatches = new ArrayList<>();
        for (int t = 0; t < n; t++) {
            int xT = xsAll.get(t);
            BigInteger yT = ysAll.get(t);
            BigInteger yHat = evalAtXFromSubset(xsBest, ysBest, xT, p);
            if (!yHat.equals(yT)) {
                BigInteger yOrig = data.shares.get(t).y;
                mismatches.add(String.format("x=%d, given_y=%s, expected_mod_p=%s",
                        xT, yOrig.toString(), yHat.toString()));
            }
        }

        System.out.println("Prime used (p): " + p);
        System.out.println("Degree assumed (k-1): " + (k - 1));
        System.out.println("Max agreement: " + best.agree + " out of " + n);
        System.out.println("Secret f(0) mod p: " + best.secret);

        if (best.ties.size() > 1) {
            System.out.println("Note: multiple subsets tied with the same agreement (" + best.ties.size() + ").");
        }

        if (mismatches.isEmpty()) {
            System.out.println("All shares are consistent with the best polynomial (no wrong shares detected).");
        } else {
            System.out.println("Suspected wrong shares (" + mismatches.size() + "):");
            for (String s : mismatches) System.out.println("  - " + s);
        }
    }
}
