class Solution:
    def lengthOfLongestSubstring(self, s: str) -> int:
        if len(s) == 0:
            return 0
        if len(s) == 1:
            return 1
        max_len = 0
        for i in range(len(s)):
            print ("i", i)
            cur_len = 1
            for j in range(i + 1, len(s)):
                print ("j",j)
                if s[j] == s[i]:
                    break
                else:
                    cur_len += 1

            if cur_len > max_len:
                max_len = cur_len
        return max_len

s=Solution()
print (s.lengthOfLongestSubstring('abcabcbb'))