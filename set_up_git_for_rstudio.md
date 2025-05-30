Setting Up Git and GitHub for R Users
================
Your Name
2025-04-15

# Git and GitHub Setup for R Users

This tutorial will guide you through setting up Git and GitHub for use
with RStudio. By following these steps, you’ll be able to version
control your R projects and collaborate with others through GitHub

## Prerequisites

Before starting, ensure you have:

- RStudio installed (latest version recommended)
- R installed (latest version recommended)
- Internet connection to download packages and access GitHub

## Installation Process

### 1. Sign up for GitHub

If you don’t already have a GitHub account, you’ll need to create one:

1.  Go to [GitHub’s join page](https://github.com/join)
2.  Follow the prompts to create your account
3.  Verify your email address

GitHub is a web-based hosting service for Git repositories, offering
collaboration features, issue tracking, and more.

### 2. Install Git

Git is a version control system that tracks changes in your code over
time. Install Git based on your operating system:

**For Windows users:**

``` r
# Download and install from:
browseURL("https://git-scm.com/download/win")
```

**For Mac users:**

``` r
# Download and install from:
browseURL("https://git-scm.com/download/mac")
```

**For Linux users:**

``` r
# Download and install from:
browseURL("https://git-scm.com/download/linux")

# Or use your distribution's package manager, e.g., for Fedora:
# sudo dnf install git-all
```

After installation, restart RStudio to ensure it recognizes Git.

### 3. Configure Git with RStudio

Now we’ll configure Git to work with RStudio by setting your user
information and authentication.

#### Set your user name and email

Use the `usethis` package to configure your Git username and email:

``` r
# Install usethis if needed
if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}

# Set your Git username and email
# Replace with your actual name and email
usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")
```

Your Git commits will be associated with this name and email.

#### Create a Personal Access Token (PAT)

GitHub requires a Personal Access Token for secure authentication:

``` r
# This will open GitHub in your browser to create a token
usethis::create_github_token()

# For older versions of usethis (< 2.0.0)
# usethis::browse_github_token()
```

Follow the prompts on GitHub to generate your token. Make sure to: 1.
Give it a descriptive name (e.g., “R-Studio-Access”) 2. Set an
appropriate expiration date 3. Select scopes (permissions) - typically
“repo” and “workflow” are sufficient 4. Copy the token immediately after
creation (you won’t be able to see it again!)

#### For Linux Users: Extend Credential Cache

Linux users may want to extend the credential cache timeout:

``` r
# Set cache timeout to approximately 30 days
usethis::use_git_config(credential.helper = "cache --timeout=2600000")
```

#### Store your Personal Access Token

There are two ways to store your PAT:

**Option 1:** Using the credentials package:

``` r
# Install credentials if needed
if (!requireNamespace("credentials", quietly = TRUE)) {
  install.packages("credentials")
}

# Replace "YourPAT" with your actual token
credentials::set_github_pat("YourPAT")
```

**Option 2:** Store in your `.Renviron` file:

``` r
# Open .Renviron file
usethis::edit_r_environ()

# In the file that opens, add this line (with your actual token):
# GITHUB_PAT=xxxyyyzzz
# Make sure the file ends with a newline
```

### 4. Restart R

After configuring Git, restart your R session:

``` r
# You can use this function:
.rs.restartR()

# Or click Session > Restart R in RStudio
```

### 5. Verify Settings

Check that everything is configured correctly:

``` r
# Check Git configuration
usethis::git_sitrep()
```

Your output should include: - Your correct username and email -
Confirmation that your PAT was found
(`Personal access token: '<found in env var>'`)

If you see any errors, carefully read the output for troubleshooting
guidance.

## Troubleshooting

If you encounter issues:

1.  **PAT not found**: Edit your `.Renviron` file manually with
    `usethis::edit_r_environ()` and ensure your token is stored
    correctly as `GITHUB_PAT=your_token_here`.

2.  **Git not detected**: Make sure Git is installed and restart
    RStudio.

3.  **Authentication failures**: Generate a new PAT if your current one
    has expired.

## Next Steps

Now that you’ve set up Git and GitHub with R:

1.  Create your first Git repository with `usethis::use_git()`
2.  Connect it to GitHub with `usethis::use_github()`
3.  Learn basic Git commands like commit, push, and pull

## References

- Gist for [Configure GitHub for
  Rstudio](https://gist.github.com/Z3tt/3dab3535007acf108391649766409421)
- Short video playlist for \[How to Use Git/GitHub with R\]
  (<https://www.youtube.com/watch?v=k-57UBYnRvw&list=PLSviU861UtD81AuyYb3SbndmAA_qTCoLe>)
- [Happy Git with R](https://happygitwithr.com/)
- [usethis package documentation](https://usethis.r-lib.org/)
- [GitHub documentation](https://docs.github.com/)

------------------------------------------------------------------------

**Congratulations!** You’ve successfully configured Git and GitHub to
work with RStudio. This setup will allow you to use version control for
your R projects and collaborate efficiently with others.
